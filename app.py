
from __future__ import annotations

import json
import io
import os
import base64
import re
import zipfile
import xml.etree.ElementTree as ET

from datetime import datetime, timezone
from typing import Callable
from urllib.parse import quote, urlparse

import pandas as pd
import requests
import streamlit as st
from Bio import Entrez


APP_TITLE = "PMC Literature Finder"
APP_SUBTITLE = "Pretraga PubMed Central baze i izvoz bibliografskih referenci"

PMC_S3_BUCKET_URL = "https://pmc-oa-opendata.s3.amazonaws.com/"

MAX_SINGLE_PDF_BYTES = 60 * 1024 * 1024
MAX_PDF_PACKAGE_BYTES = 200 * 1024 * 1024

REQUEST_TIMEOUT = (10, 90)


# ============================================================
# 1. PMC IDENTIFIERS AND METADATA
# ============================================================

def normalize_pmc_id(value: str) -> str:
    value = str(value or "").strip()

    if value and not value.upper().startswith("PMC"):
        return f"PMC{value}"

    return value


def extract_summary_records(
    summary_json: dict,
    requested_ids: list[str]
) -> list[dict[str, str]]:
    """Normalize PMC ESummary JSON metadata."""

    result = summary_json.get("result", {})

    returned_uids = [
        str(uid) for uid in result.get("uids", [])
    ]

    ordered_uids = returned_uids or [
        str(uid) for uid in requested_ids
    ]

    records = []

    for uid in ordered_uids:

        item = result.get(uid, {})

        if not item:
            continue

        article_ids = (
            item.get("articleids", [])
            or item.get("article_ids", [])
            or []
        )

        pmcid = ""
        pmid = ""
        doi = ""

        for article_id in article_ids:

            if not isinstance(article_id, dict):
                continue

            id_type = str(
                article_id.get("idtype", "")
            ).lower()

            value = str(
                article_id.get("value", "")
            ).strip()

            if id_type in {"pmc", "pmcid"}:
                pmcid = normalize_pmc_id(value)

            elif id_type in {"pubmed", "pmid"}:
                pmid = value

            elif id_type == "doi":
                doi = value

        pmcid = pmcid or normalize_pmc_id(uid)

        if not doi:
            doi = str(item.get("doi", "") or "")

        authors_data = (
            item.get("authors", [])
            or item.get("authorlist", [])
            or item.get("AuthorList", [])
            or []
        )

        authors = []

        for author in authors_data:

            if isinstance(author, dict):
                name = str(
                    author.get(
                        "name",
                        author.get("Name", "")
                    )
                ).strip()

            else:
                name = str(author).strip()

            if name:
                authors.append(name)

        pubdate = str(
            item.get(
                "pubdate",
                item.get("PubDate", "")
            ) or ""
        ).strip()

        year_match = re.search(
            r"\b(?:19|20)\d{2}\b",
            pubdate
        )

        year = (
            year_match.group(0)
            if year_match
            else ""
        )

        pages_raw = str(
            item.get(
                "pages",
                item.get("Pages", "")
            ) or ""
        ).strip()

        first_page = ""
        last_page = ""

        if pages_raw:

            page_match = re.match(
                r"^([A-Za-z0-9]+)\s*[-–]\s*([A-Za-z0-9]+)$",
                pages_raw
            )

            if page_match:
                first_page, last_page = page_match.groups()

            else:
                first_page = pages_raw

        pmc_url = (
            f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
            if pmcid
            else ""
        )

        records.append({

            "PMCID": pmcid,
            "PMID": pmid,

            "Title": str(
                item.get("title", item.get("Title", ""))
                or ""
            ).strip(),

            "Authors": "; ".join(authors),

            "Year": year,
            "Publication date": pubdate,

            "Journal": str(
                item.get(
                    "fulljournalname",
                    item.get("source", item.get("Source", ""))
                ) or ""
            ).strip(),

            "Volume": str(
                item.get("volume", item.get("Volume", ""))
                or ""
            ).strip(),

            "Issue": str(
                item.get("issue", item.get("Issue", ""))
                or ""
            ).strip(),

            "First page": first_page,
            "Last page": last_page,

            "Article locator": str(
                item.get(
                    "elocationid",
                    item.get("elocationId", "")
                ) or ""
            ).strip(),

            "DOI": doi,
            "PMC URL": pmc_url,

        })

    return records


# ============================================================
# 2. CITATION FORMATTING
# ============================================================

def apa_author_name(author: str) -> str:

    author = author.strip()

    if not author:
        return author

    if "," in author:
        surname, given = [
            part.strip()
            for part in author.split(",", 1)
        ]

    else:
        parts = author.split()

        if len(parts) < 2:
            return author

        surname = " ".join(parts[:-1])
        given = parts[-1]

    initials = []

    for token in re.split(r"[\s.]+", given):

        if not token:
            continue

        chunks = token.split("-")
        chunk_initials = []

        for chunk in chunks:

            if not chunk:
                continue

            letters = (
                list(chunk)
                if chunk.isalpha()
                and chunk.isupper()
                and len(chunk) > 1
                else [chunk[0]]
            )

            chunk_initials.append(
                " ".join(
                    f"{letter.upper()}."
                    for letter in letters
                )
            )

        if chunk_initials:
            initials.append("-".join(chunk_initials))

    return (
        f"{surname}, {' '.join(initials)}"
        if initials
        else author
    )


def apa_authors(authors: list[str]) -> str:

    if not authors:
        return ""

    formatted = [
        apa_author_name(author)
        for author in authors
    ]

    if len(formatted) <= 20:
        shown = formatted

    else:
        shown = (
            formatted[:19]
            + ["…"]
            + [formatted[-1]]
        )

    if len(shown) == 1:
        return shown[0]

    if len(shown) == 2:
        return f"{shown[0]} & {shown[1]}"

    return (
        ", ".join(shown[:-1])
        + f", & {shown[-1]}"
    )


def vancouver_authors(authors: list[str]) -> str:

    if not authors:
        return ""

    shown = authors[:6]

    result = ", ".join(shown)

    if len(authors) > 6:
        result += ", et al."

    return result


def format_reference(
    record: dict[str, str],
    style: str
) -> str:

    authors = [
        author.strip()
        for author in record.get("Authors", "").split(";")
        if author.strip()
    ]

    year = record.get("Year") or "n.d."
    title = record.get("Title") or "[Untitled article]"

    journal = record.get("Journal", "")
    volume = record.get("Volume", "")
    issue = record.get("Issue", "")

    first_page = record.get("First page", "")
    last_page = record.get("Last page", "")

    locator = record.get("Article locator", "")
    doi = record.get("DOI", "")

    doi_text = (
        f"https://doi.org/{doi}"
        if doi
        else ""
    )

    pages = ""

    if first_page and last_page:

        pages = (
            first_page
            if first_page == last_page
            else f"{first_page}–{last_page}"
        )

    elif first_page:
        pages = first_page

    elif locator:
        pages = locator

    volume_issue = volume

    if issue:
        volume_issue += (
            f"({issue})"
            if volume
            else f"({issue})"
        )

    if style == "APA 7":

        author_text = (
            apa_authors(authors)
            or "Unknown author"
        )

        journal_part = journal

        if volume_issue:
            journal_part += f", {volume_issue}"

        if pages:
            journal_part += f", {pages}"

        if journal_part:

            result = (
                f"{author_text} ({year}). "
                f"{title}. {journal_part}."
            )

        else:

            result = (
                f"{author_text} ({year}). {title}."
            )

        if doi_text:
            result += f" {doi_text}"

        return result

    if style == "Vancouver":

        author_text = vancouver_authors(authors)

        result = (
            f"{author_text}. {title}. "
            f"{journal}. {year}"
        )

        if volume:

            result += f";{volume}"

            if issue:
                result += f"({issue})"

        if pages:
            result += f":{pages}"

        result += "."

        if doi:
            result += f" doi:{doi}."

        return result

    # Harvard

    author_text = (
        apa_authors(authors)
        or "Unknown author"
    )

    result = (
        f"{author_text} ({year}) "
        f"'{title}', {journal}"
    )

    if volume_issue:
        result += f", vol. {volume_issue}"

    if pages:
        result += f", pp. {pages}"

    result += "."

    if doi_text:
        result += f" Available at: {doi_text}."

    return result


# ============================================================
# 3. BIBTEX AND RIS EXPORT
# ============================================================

def bibtex_escape(value: str) -> str:

    return (
        value.replace("\\", r"\textbackslash ")
        .replace("{", r"\{")
        .replace("}", r"\}")
    )


def make_bibtex(
    records: list[dict[str, str]]
) -> str:

    entries = []

    for i, record in enumerate(records, start=1):

        authors = [
            author.strip()
            for author in record.get("Authors", "").split(";")
            if author.strip()
        ]

        first_surname = (
            authors[0].split(",")[0].split()[-1]
            if authors
            else "article"
        )

        year = record.get("Year") or "nd"
        suffix = record.get("PMCID") or str(i)

        key = re.sub(
            r"[^A-Za-z0-9]+",
            "",
            f"{first_surname}{year}{suffix}"
        )

        fields = {

            "title": record.get("Title", ""),

            "author": " and ".join(authors),

            "journal": record.get("Journal", ""),

            "year": record.get("Year", ""),

            "volume": record.get("Volume", ""),

            "number": record.get("Issue", ""),

            "pages": (
                f"{record.get('First page', '')}--"
                f"{record.get('Last page', '')}"
                if record.get("First page")
                and record.get("Last page")
                else record.get("Article locator", "")
            ),

            "doi": record.get("DOI", ""),

            "url": record.get("PMC URL", ""),

        }

        if authors:

            fields["author"] = " and ".join(
                apa_author_name(author).replace(". ", " ")
                for author in authors
            )

        body = [
            f"  {name} = {{{bibtex_escape(value)}}}"
            for name, value in fields.items()
            if value
        ]

        entries.append(
            "@article{"
            + key
            + ",\n"
            + ",\n".join(body)
            + "\n}"
        )

    return "\n\n".join(entries)


def make_ris(
    records: list[dict[str, str]]
) -> str:

    lines = []

    for record in records:

        lines.append("TY  - JOUR")

        if record.get("Title"):
            lines.append(f"TI  - {record['Title']}")

        authors = [
            author.strip()
            for author in record.get("Authors", "").split(";")
            if author.strip()
        ]

        for author in authors:
            lines.append(f"AU  - {author}")

        fields = [

            ("Journal", "JO"),
            ("Year", "PY"),
            ("Volume", "VL"),
            ("Issue", "IS"),
            ("First page", "SP"),
            ("Last page", "EP"),
            ("DOI", "DO"),
            ("PMC URL", "UR"),

        ]

        for field, tag in fields:

            if record.get(field):
                lines.append(
                    f"{tag}  - {record[field]}"
                )

        if record.get("PMCID"):
            lines.append(
                f"AN  - {record['PMCID']}"
            )

        lines.append("ER  - ")
        lines.append("")

    return "\n".join(lines)


# ============================================================
# 4. ZIP PACKAGE
# ============================================================

def create_zip(
    records: list[dict[str, str]],
    selected_style: str,
    query: str,
    total_count: str,
    pdf_files: list[dict[str, object]] | None = None,
    pdf_report: list[dict[str, str]] | None = None,
) -> bytes:

    pdf_files = pdf_files or []
    pdf_report = pdf_report or []

    df = pd.DataFrame(records)

    citations_selected = "\n\n".join(
        format_reference(record, selected_style)
        for record in records
    )

    all_styles = {

        style: "\n\n".join(
            format_reference(record, style)
            for record in records
        )

        for style in [
            "APA 7",
            "Vancouver",
            "Harvard"
        ]

    }

    buffer = io.BytesIO()

    with zipfile.ZipFile(
        buffer,
        "w",
        compression=zipfile.ZIP_DEFLATED
    ) as archive:

        archive.writestr(
            "references.csv",
            df.to_csv(index=False).encode("utf-8-sig")
        )

        archive.writestr(
            "references.ris",
            make_ris(records).encode("utf-8")
        )

        archive.writestr(
            "references.bib",
            make_bibtex(records).encode("utf-8")
        )

        archive.writestr(
            "citations_selected_style.txt",
            citations_selected.encode("utf-8")
        )

        for style, text in all_styles.items():

            filename = (
                f"citations_{style.lower().replace(' ', '_')}.txt"
            )

            archive.writestr(
                filename,
                text.encode("utf-8")
            )

        for pdf in pdf_files:

            filename = str(
                pdf.get("filename", "article.pdf")
            )

            data = pdf.get("data", b"")

            if isinstance(data, bytes) and data:

                archive.writestr(
                    f"pdfs/{filename}",
                    data
                )

        if pdf_report:

            report_df = pd.DataFrame(pdf_report)

            archive.writestr(
                "pdf_download_report.csv",
                report_df.to_csv(
                    index=False
                ).encode("utf-8-sig")
            )

        pdf_note = (

            f"PDF files included: {len(pdf_files)}\n"

            "PDF retrieval was limited to the PMC Open Access Subset "
            "and records for which PMC reports a PDF file. "

            "The license for each downloaded article is listed in "
            "pdf_download_report.csv. "

            "Downloaded PDFs remain subject to their individual "
            "license terms; inclusion in this ZIP does not grant "
            "additional reuse rights.\n"

            if pdf_report

            else "PDF files included: 0 (bibliographic package only)\n"

        )

        info = (

            f"Source: PubMed Central (PMC)\n"

            f"Query: {query}\n"

            f"Total matching records reported by PMC: {total_count}\n"

            f"Records retrieved: {len(records)}\n"

            f"PDF files included: {len(pdf_files)}\n"

            f"Export date (UTC): "
            f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}\n"

            "Note: citation strings are automatically formatted "
            "from metadata returned by PMC; verify each reference "
            "against the journal's style requirements.\n"

            + pdf_note

        )

        archive.writestr(
            "search_info.txt",
            info.encode("utf-8")
        )

    buffer.seek(0)

    return buffer.getvalue()


# ============================================================
# 5. PMC SEARCH
# ============================================================

def search_pmc(
    query: str,
    max_results: int,
    email: str,
    api_key: str = ""
) -> tuple[list[dict[str, str]], str]:

    Entrez.email = email.strip()
    Entrez.tool = "PMC_Literature_Finder_Streamlit"
    Entrez.api_key = api_key.strip() or None

    with Entrez.esearch(
        db="pmc",
        term=query,
        retmax=max_results,
        retmode="xml"
    ) as handle:

        search_result = Entrez.read(handle)

    ids = [
        str(item)
        for item in search_result.get("IdList", [])
    ]

    total_count = str(
        search_result.get("Count", "0")
    )

    if not ids:
        return [], total_count

    with Entrez.esummary(
        db="pmc",
        id=",".join(ids),
        retmode="json"
    ) as handle:

        summary_json = json.load(handle)

    records = extract_summary_records(
        summary_json,
        ids
    )

    return records, total_count


# ============================================================
# 6. PMC OPEN ACCESS PDF SEARCH
# ============================================================

def search_pmc_pdf_ids(
    query: str,
    max_results: int,
    email: str,
    api_key: str = "",
) -> tuple[set[str], str]:

    """Return PMC IDs matching OA-subset and has_pdf filters."""

    Entrez.email = email.strip()
    Entrez.tool = "PMC_Literature_Finder_Streamlit"
    Entrez.api_key = api_key.strip() or None

    pdf_query = (
        f"({query}) AND open_access[filter] AND has_pdf[filter]"
    )

    with Entrez.esearch(
        db="pmc",
        term=pdf_query,
        retmax=max_results,
        retmode="xml"
    ) as handle:

        result = Entrez.read(handle)

    ids = {
        str(value)
        for value in result.get("IdList", [])
    }

    return ids, str(result.get("Count", "0"))


# ============================================================
# 7. PMC CLOUD PDF METADATA
# ============================================================

def s3_uri_to_https(value: str) -> str:

    """Convert public PMC S3 URI to HTTPS."""

    parsed = urlparse(str(value or ""))

    if parsed.scheme == "https":
        return str(value)

    if (
        parsed.scheme != "s3"
        or parsed.netloc != "pmc-oa-opendata"
    ):

        raise ValueError(
            "PMC je u metapodacima vratio neočekivanu PDF adresu."
        )

    key = quote(
        parsed.path.lstrip("/"),
        safe="/"
    )

    url = (
        f"https://{parsed.netloc}.s3.amazonaws.com/{key}"
    )

    return (
        f"{url}?{parsed.query}"
        if parsed.query
        else url
    )


def find_pmc_pdf_metadata(
    pmcid: str,
    session: requests.Session
) -> dict[str, object] | None:

    """Find newest PMC Cloud JSON metadata object with OA PDF URL."""

    pmcid = normalize_pmc_id(pmcid)

    params = {

        "list-type": "2",
        "prefix": f"metadata/{pmcid}.",
        "max-keys": "100",

    }

    response = session.get(
        PMC_S3_BUCKET_URL,
        params=params,
        timeout=REQUEST_TIMEOUT
    )

    response.raise_for_status()

    root = ET.fromstring(response.content)

    keys = [

        (element.text or "").strip()

        for element in root.iter()

        if element.tag.rsplit("}", 1)[-1] == "Key"

    ]

    metadata_records = []

    for key in keys:

        if not key.endswith(".json"):
            continue

        metadata_url = (
            PMC_S3_BUCKET_URL
            + quote(key, safe="/")
        )

        metadata_response = session.get(
            metadata_url,
            timeout=REQUEST_TIMEOUT
        )

        metadata_response.raise_for_status()

        item = metadata_response.json()

        if (
            str(item.get("pmcid", "")).upper()
            != pmcid.upper()
        ):
            continue

        if item.get("is_pmc_openaccess") is not True:
            continue

        if not item.get("pdf_url"):
            continue

        metadata_records.append(item)

    if not metadata_records:
        return None

    return max(
        metadata_records,
        key=lambda item: int(
            item.get("version", 0) or 0
        )
    )


# ============================================================
# 8. DOWNLOAD PMC PDF FILES
# ============================================================

def download_pmc_pdfs(
    records: list[dict[str, str]],
    query: str,
    email: str,
    api_key: str = "",
    progress_callback: Callable[
        [int, int, str], None
    ] | None = None,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, str]],
    str
]:

    """Download PDFs for search hits in PMC OA subset."""

    Entrez.email = email.strip()
    Entrez.tool = "PMC_Literature_Finder_Streamlit"
    Entrez.api_key = api_key.strip() or None

    selected_numeric_ids = {

        normalize_pmc_id(
            record.get("PMCID", "")
        ).replace("PMC", "", 1)

        for record in records

    }

    eligible_ids, _eligible_total = search_pmc_pdf_ids(
        query,
        max_results=10000,
        email=email,
        api_key=api_key
    )

    eligible_ids.intersection_update(
        selected_numeric_ids
    )

    eligible_count = str(
        len(eligible_ids)
    )

    pdf_files = []
    report = []

    downloaded_bytes = 0

    total_steps = max(
        len(records),
        1
    )

    with requests.Session() as session:

        session.headers.update({

            "User-Agent": (
                "PMC-Literature-Finder/1.0 "
                "(academic literature retrieval)"
            )

        })

        for index, record in enumerate(
            records,
            start=1
        ):

            pmcid = normalize_pmc_id(
                record.get("PMCID", "")
            )

            numeric_id = pmcid.replace(
                "PMC",
                "",
                1
            )

            if progress_callback:

                progress_callback(
                    index,
                    total_steps,
                    pmcid
                )

            if numeric_id not in eligible_ids:

                report.append({

                    "PMCID": pmcid,

                    "Title": record.get("Title", ""),

                    "Status": "Nije u PMC OA + has_pdf skupu",

                    "License": "",

                    "PDF file": "",

                    "Details": (
                        "PMC pretraga nije označila ovaj zapis "
                        "kao OA članak sa dostupnim PDF-om."
                    ),

                })

                continue

            try:

                metadata = find_pmc_pdf_metadata(
                    pmcid,
                    session
                )

                if not metadata:

                    report.append({

                        "PMCID": pmcid,

                        "Title": record.get("Title", ""),

                        "Status": "PDF metapodaci nisu pronađeni",

                        "License": "",

                        "PDF file": "",

                        "Details": (
                            "Zapis je prošao PMC filter, ali "
                            "odgovarajući PMC Cloud JSON/PDF URL "
                            "nije pronađen."
                        ),

                    })

                    continue

                license_code = str(
                    metadata.get("license_code")
                    or "Nije navedena u mašinski čitljivom polju"
                )

                pdf_url = s3_uri_to_https(
                    str(metadata.get("pdf_url", ""))
                )

                response = session.get(
                    pdf_url,
                    stream=True,
                    timeout=REQUEST_TIMEOUT
                )

                response.raise_for_status()

                content_length = response.headers.get(
                    "Content-Length"
                )

                if (
                    content_length
                    and int(content_length) > MAX_SINGLE_PDF_BYTES
                ):

                    response.close()

                    report.append({

                        "PMCID": pmcid,

                        "Title": record.get("Title", ""),

                        "Status": "Preskočeno: PDF je prevelik",

                        "License": license_code,

                        "PDF file": "",

                        "Details": (
                            "Veličina prelazi limit od "
                            f"{MAX_SINGLE_PDF_BYTES // (1024 * 1024)} MB po radu."
                        ),

                    })

                    continue

                chunks = []
                current_size = 0
                too_large = False

                for chunk in response.iter_content(
                    chunk_size=1024 * 256
                ):

                    if not chunk:
                        continue

                    current_size += len(chunk)

                    if current_size > MAX_SINGLE_PDF_BYTES:
                        too_large = True
                        break

                    chunks.append(chunk)

                response.close()

                if too_large:

                    report.append({

                        "PMCID": pmcid,

                        "Title": record.get("Title", ""),

                        "Status": "Preskočeno: PDF je prevelik",

                        "License": license_code,

                        "PDF file": "",

                        "Details": (
                            "Veličina prelazi limit od "
                            f"{MAX_SINGLE_PDF_BYTES // (1024 * 1024)} MB po radu."
                        ),

                    })

                    continue

                pdf_data = b"".join(chunks)

                if not pdf_data.startswith(b"%PDF"):

                    raise ValueError(
                        "Preuzeti sadržaj nije prepoznat kao PDF dokument."
                    )

                if (
                    downloaded_bytes + len(pdf_data)
                    > MAX_PDF_PACKAGE_BYTES
                ):

                    report.append({

                        "PMCID": pmcid,

                        "Title": record.get("Title", ""),

                        "Status": (
                            "Preskočeno: dostignut limit ZIP paketa"
                        ),

                        "License": license_code,

                        "PDF file": "",

                        "Details": (
                            "Limit ukupnog paketa je "
                            f"{MAX_PDF_PACKAGE_BYTES // (1024 * 1024)} MB."
                        ),

                    })

                    continue

                version = str(
                    metadata.get("version", "1")
                )

                filename = (
                    f"{pmcid}.v{version}.pdf"
                )

                pdf_files.append({

                    "filename": filename,

                    "data": pdf_data,

                })

                downloaded_bytes += len(pdf_data)

                report.append({

                    "PMCID": pmcid,

                    "Title": record.get("Title", ""),

                    "Status": "Preuzet",

                    "License": license_code,

                    "PDF file": filename,

                    "Details": (
                        "PDF preuzet iz zvaničnog PMC Cloud skupa."
                    ),

                })

            except Exception as exc:

                report.append({

                    "PMCID": pmcid,

                    "Title": record.get("Title", ""),

                    "Status": "Greška pri preuzimanju",

                    "License": "",

                    "PDF file": "",

                    "Details": (
                        f"{type(exc).__name__}: {exc}"
                    ),

                })

    return pdf_files, report, eligible_count


# ============================================================
# 9. APPLICATION BACKGROUND AND STYLING
# ============================================================

def set_background(
    image_file: str = "background.jpg"
) -> None:

    """Use PCR Primer Checker background when supplied."""

    if os.path.exists(image_file):

        with open(image_file, "rb") as image:

            encoded = base64.b64encode(
                image.read()
            ).decode("utf-8")

        background_css = (
            f'url("data:image/jpeg;base64,{encoded}")'
        )

        size_css = "cover"

    else:

        background_css = (
            "linear-gradient(135deg, "
            "#152b45 0%, #315c72 48%, #78939a 100%)"
        )

        size_css = "cover"

    st.markdown(
        f"""
        <style>

        .stApp {{

            background-image: {background_css};
            background-size: {size_css};
            background-position: center;
            background-attachment: fixed;

        }}

        [data-testid="stMainBlockContainer"] {{

            background-color: rgba(255, 255, 255, 0.88);

            padding: 2rem 2.2rem;

            border-radius: 20px;

            box-shadow: 0 0 30px rgba(0, 0, 0, 0.28);

            max-width: 1250px;

        }}

        [data-testid="stMainBlockContainer"] h1,
        [data-testid="stMainBlockContainer"] h2,
        [data-testid="stMainBlockContainer"] h3,
        [data-testid="stMainBlockContainer"] p,
        [data-testid="stMainBlockContainer"] label,
        [data-testid="stMainBlockContainer"] [data-testid="stCaptionContainer"] {{

            color: #17212b;

        }}

        [data-testid="stMainBlockContainer"] [data-testid="stMetric"] {{

            background: rgba(235, 242, 246, 0.92);

            border: 1px solid rgba(35, 65, 85, 0.12);

            padding: 0.8rem 1rem;

            border-radius: 12px;

        }}

        [data-testid="stMainBlockContainer"] [data-testid="stDataFrame"] {{

            border: 1px solid rgba(35, 65, 85, 0.16);

            border-radius: 10px;

            overflow: hidden;

        }}

        [data-testid="stMainBlockContainer"] .stButton > button,
        [data-testid="stMainBlockContainer"] .stDownloadButton > button {{

            border-radius: 9px;

            font-weight: 600;

        }}

        [data-testid="stMainBlockContainer"] [data-testid="stForm"] {{

            background: rgba(246, 249, 251, 0.88);

            border: 1px solid rgba(35, 65, 85, 0.13);

            border-radius: 14px;

            padding: 1.2rem 1.2rem 0.8rem 1.2rem;

        }}

        @media (prefers-color-scheme: dark) {{

            [data-testid="stMainBlockContainer"] {{

                background-color: rgba(27, 34, 42, 0.91);

            }}

            [data-testid="stMainBlockContainer"] h1,
            [data-testid="stMainBlockContainer"] h2,
            [data-testid="stMainBlockContainer"] h3,
            [data-testid="stMainBlockContainer"] p,
            [data-testid="stMainBlockContainer"] label,
            [data-testid="stMainBlockContainer"] [data-testid="stCaptionContainer"] {{

                color: #f0f3f6;

            }}

            [data-testid="stMainBlockContainer"] [data-testid="stMetric"],
            [data-testid="stMainBlockContainer"] [data-testid="stForm"] {{

                background: rgba(45, 55, 66, 0.94);

                border-color: rgba(220, 230, 240, 0.15);

            }}

        }}

        </style>
        """,
        unsafe_allow_html=True,
    )


# ============================================================
# 10. MAIN STREAMLIT APPLICATION
# ============================================================

def main() -> None:

    st.set_page_config(
        page_title=APP_TITLE,
        page_icon="📚",
        layout="wide"
    )

    set_background()

    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    st.markdown(
        """
        <div style="text-align:center; padding: 0.2rem 0 1rem 0;">

          <h1 style="margin-bottom:0.25rem;">
            📚 PMC Literature Finder
          </h1>

          <h3 style="margin-top:0; font-weight:500;">
            Pretraga i organizacija naučne literature
          </h3>

          <p style="margin-bottom:0.15rem;">
            Laboratorija za Molekularnu biologiju
          </p>

          <p style="margin:0.1rem;">
            Departman za biologiju i ekologiju ·
            Prirodno-matematički fakultet
          </p>

          <p style="margin:0.1rem;">
            Univerzitet u Nišu
          </p>

        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div style="
            height:3px;
            background:linear-gradient(
                90deg,
                #4c7187,
                #b4cbd2,
                #4c7187
            );
            border-radius:3px;
            margin:0.4rem 0 1.2rem 0;
        "></div>
        """,
        unsafe_allow_html=True,
    )

    # --------------------------------------------------------
    # Application information
    # --------------------------------------------------------

    with st.expander(
        "O aplikaciji i načinu pretrage",
        expanded=False
    ):

        st.markdown(
            "Aplikacija pretražuje **PubMed Central (PMC)** "
            "preko NCBI E-utilities, preuzima bibliografske "
            "metapodatke i priprema CSV, RIS, BibTeX i "
            "formatirane citate. "
            "Po izboru, PDF-ovi se preuzimaju iz zvaničnog "
            "PMC Cloud skupa samo za radove iz Open Access "
            "Subset-a za koje PMC prijavljuje dostupnu PDF datoteku."
        )

        st.markdown(
            "**Primeri upita:** "
            "`glioblastoma AND IDH1` · "
            "`\"honey bee\" AND microbiome` · "
            "`(cancer OR glioma) AND methylation`"
        )

        st.caption(
            "Obavezna je e-mail adresa za NCBI E-utilities. "
            "API ključ je opcioni. "
            "PDF preuzimanje je dodatni, odvojeni korak "
            "i zavisi od dostupnosti PDF-a i licence konkretnog rada."
        )

    # --------------------------------------------------------
    # Search form
    # --------------------------------------------------------

    st.subheader("Pretraga literature")

    with st.form(
        "pmc_search_form",
        clear_on_submit=False
    ):

        query = st.text_input(
            "Ključne reči / PMC upit",
            value="glioblastoma AND IDH1",
            help=(
                "Možeš koristiti AND, OR, NOT, navodnike "
                "i PMC/Entrez oznake polja."
            ),
        )

        col1, col2 = st.columns(2)

        with col1:

            max_results = st.slider(
                "Broj rezultata",
                min_value=1,
                max_value=100,
                value=10,
                step=1
            )

            citation_style = st.selectbox(
                "Stil citiranja",
                [
                    "APA 7",
                    "Vancouver",
                    "Harvard"
                ]
            )

        with col2:

            email = st.text_input(
                "E-mail za NCBI E-utilities",
                placeholder="ime.prezime@institucija.rs"
            )

            api_key = st.text_input(
                "NCBI API key (opciono)",
                type="password"
            )

        submitted = st.form_submit_button(
            "🔎 Pretraži PMC",
            type="primary",
            use_container_width=True
        )

    # --------------------------------------------------------
    # Execute search
    # --------------------------------------------------------

    if submitted:

        if not query.strip():

            st.error(
                "Unesi ključne reči ili PMC upit."
            )

            return

        if not email.strip() or "@" not in email:

            st.error(
                "Unesi važeću e-mail adresu. "
                "NCBI je traži za identifikaciju automatizovanih zahteva."
            )

            return

        with st.spinner(
            "Pretražujem PMC i preuzimam bibliografske podatke..."
        ):

            try:

                records, total_count = search_pmc(
                    query.strip(),
                    max_results,
                    email,
                    api_key
                )

            except Exception as exc:

                st.error(
                    f"Pretraga nije uspela "
                    f"({type(exc).__name__}): {exc}"
                )

                st.info(
                    "Proveri internet konekciju, e-mail adresu "
                    "i sintaksu upita. Pokušaj i sa jednostavnijim ključnim rečima."
                )

                return

        if not records:

            st.warning(
                f"PMC nije vratio zapise za upit: {query}"
            )

            st.session_state.pop(
                "pmc_records",
                None
            )

            return

        st.session_state["pmc_records"] = records

        st.session_state["pmc_query"] = query.strip()

        st.session_state["pmc_total_count"] = total_count

        st.session_state["pmc_style"] = citation_style

        st.session_state["pmc_email"] = email.strip()

        st.session_state["pmc_api_key"] = api_key.strip()

        # Clear previous PDF results after a new search.

        st.session_state.pop(
            "pmc_pdf_files",
            None
        )

        st.session_state.pop(
            "pmc_pdf_report",
            None
        )

        st.session_state.pop(
            "pmc_pdf_eligible_count",
            None
        )

    # --------------------------------------------------------
    # Display search results
    # --------------------------------------------------------

    records = st.session_state.get(
        "pmc_records"
    )

    if not records:

        st.info(
            "Unesi ključne reči i e-mail adresu, "
            "pa pokreni pretragu."
        )

        return

    st.markdown("---")

    st.subheader("Rezultati pretrage")

    metric1, metric2, metric3 = st.columns(3)

    metric1.metric(
        "PMC pogodaka",
        st.session_state.get(
            "pmc_total_count",
            "—"
        )
    )

    metric2.metric(
        "Preuzetih zapisa",
        len(records)
    )

    metric3.metric(
        "Izabrani stil",
        st.session_state.get(
            "pmc_style",
            "APA 7"
        )
    )

    st.caption(
        f"Upit: `{st.session_state.get('pmc_query', '')}`"
    )

    display_columns = [

        "PMCID",
        "PMID",
        "Title",
        "Authors",
        "Year",
        "Journal",
        "DOI"

    ]

    results_df = pd.DataFrame(
        records
    )[display_columns]

    st.dataframe(
        results_df,
        use_container_width=True,
        hide_index=True
    )

    # --------------------------------------------------------
    # Bibliographic citations
    # --------------------------------------------------------

    st.subheader("Bibliografski citati")

    selected_style = st.selectbox(
        "Format citiranja",

        [
            "APA 7",
            "Vancouver",
            "Harvard"
        ],

        index=[
            "APA 7",
            "Vancouver",
            "Harvard"
        ].index(
            st.session_state.get(
                "pmc_style",
                "APA 7"
            )
        ),

        key="display_style",
    )

    citation_text = "\n\n".join(

        format_reference(
            record,
            selected_style
        )

        for record in records

    )

    st.text_area(
        "Generisane reference",
        value=citation_text,
        height=250
    )

    # --------------------------------------------------------
    # Basic download options
    # --------------------------------------------------------

    st.subheader("Preuzimanje rezultata")

    csv_bytes = (
        pd.DataFrame(records)
        .to_csv(index=False)
        .encode("utf-8-sig")
    )

    zip_bytes = create_zip(

        records,

        selected_style,

        st.session_state.get(
            "pmc_query",
            ""
        ),

        st.session_state.get(
            "pmc_total_count",
            ""
        ),

    )

    col_csv, col_zip, col_citations = st.columns(3)

    with col_csv:

        st.download_button(

            "📊 Preuzmi CSV",

            data=csv_bytes,

            file_name="pmc_references.csv",

            mime="text/csv",

            use_container_width=True,

        )

    with col_zip:

        st.download_button(

            "📦 Preuzmi ZIP paket",

            data=zip_bytes,

            file_name="pmc_reference_package.zip",

            mime="application/zip",

            type="primary",

            use_container_width=True,

        )

    with col_citations:

        st.download_button(

            "📄 Preuzmi citate",

            data=citation_text.encode("utf-8"),

            file_name=(
                "pmc_citations_"
                f"{selected_style.lower().replace(' ', '_')}.txt"
            ),

            mime="text/plain",

            use_container_width=True,

        )

    st.caption(
        "Osnovni ZIP sadrži CSV, RIS, BibTeX, citate u APA 7/"
        "Vancouver/Harvard stilu i informacije o pretrazi. "
        "Citati se automatski formiraju iz metapodataka PMC-a; "
        "preporučuje se provera pre predaje rukopisa."
    )

    # ========================================================
    # 11. PDF DOWNLOAD SECTION
    # ========================================================

    st.markdown("---")

    st.subheader("Preuzimanje PDF radova")

    st.markdown(
        "Dodatno možeš preuzeti PDF verzije radova koji pripadaju "
        "**PMC Open Access Subset-u** i za koje PMC označava da PDF postoji. "
        "Preuzimanje se obavlja preko zvaničnog PMC Cloud servisa."
    )

    st.info(
        "Dostupnost PDF-a nije zagarantovana za svaki rezultat. "
        "ZIP će sadržati samo uspešno preuzete PDF-ove, dok će CSV "
        "izveštaj navesti radove koji nisu bili dostupni ili čije "
        "preuzimanje nije uspelo. Svaki PDF ostaje pod uslovima "
        "licence navedenim u izveštaju."
    )

    # --------------------------------------------------------
    # PDF download button
    # --------------------------------------------------------

    if st.button(

        "📥 Preuzmi dostupne PDF-ove i pripremi ZIP",

        type="secondary",

        use_container_width=True,

    ):

        progress = st.progress(0)

        status = st.empty()

        def update_progress(
            current: int,
            total: int,
            pmcid: str
        ) -> None:

            progress.progress(
                min(
                    current / max(total, 1),
                    1.0
                )
            )

            status.caption(
                f"Provera i preuzimanje: "
                f"{current}/{total} — {pmcid}"
            )

        with st.spinner(
            "Proveravam PMC Open Access dostupnost "
            "i preuzimam PDF datoteke..."
        ):

            try:

                pdf_files, pdf_report, eligible_count = (
                    download_pmc_pdfs(

                        records=records,

                        query=st.session_state.get(
                            "pmc_query",
                            ""
                        ),

                        email=st.session_state.get(
                            "pmc_email",
                            ""
                        ),

                        api_key=st.session_state.get(
                            "pmc_api_key",
                            ""
                        ),

                        progress_callback=update_progress,

                    )
                )

                st.session_state["pmc_pdf_files"] = pdf_files

                st.session_state["pmc_pdf_report"] = pdf_report

                st.session_state["pmc_pdf_eligible_count"] = (
                    eligible_count
                )

            except Exception as exc:

                st.error(
                    "Provera PDF dostupnosti nije uspela "
                    f"({type(exc).__name__}): {exc}"
                )

                st.info(
                    "Bibliografski CSV i osnovni ZIP paket "
                    "i dalje su dostupni iznad."
                )

                return

        progress.progress(1.0)

        status.success(
            "Provera PDF dostupnosti je završena."
        )

    # --------------------------------------------------------
    # PDF download report and extended ZIP
    # --------------------------------------------------------

    if "pmc_pdf_report" in st.session_state:

        pdf_files = st.session_state.get(
            "pmc_pdf_files",
            []
        )

        pdf_report = st.session_state.get(
            "pmc_pdf_report",
            []
        )

        downloaded_count = len(pdf_files)

        eligible_count = st.session_state.get(
            "pmc_pdf_eligible_count",
            "0"
        )

        pdf_col1, pdf_col2, pdf_col3 = st.columns(3)

        pdf_col1.metric(
            "PDF kandidati prema PMC filteru",
            str(eligible_count)
        )

        pdf_col2.metric(
            "Uspešno preuzetih PDF-ova",
            downloaded_count
        )

        total_pdf_size = sum(
            len(item["data"])
            for item in pdf_files
        )

        pdf_col3.metric(
            "Ukupna veličina PDF-ova",
            f"{total_pdf_size / (1024 * 1024):.1f} MB"
        )

        report_df = pd.DataFrame(
            pdf_report
        )

        st.dataframe(
            report_df,
            use_container_width=True,
            hide_index=True
        )

        # Create extended ZIP with PDF files and references.

        extended_zip = create_zip(

            records,

            selected_style,

            st.session_state.get(
                "pmc_query",
                ""
            ),

            st.session_state.get(
                "pmc_total_count",
                ""
            ),

            pdf_files=pdf_files,

            pdf_report=pdf_report,

        )

        st.download_button(

            "📦 Preuzmi kompletan ZIP sa dostupnim PDF-ovima",

            data=extended_zip,

            file_name="pmc_references_and_pdfs.zip",

            mime="application/zip",

            type="primary",

            use_container_width=True,

        )


# ============================================================
# 12. RUN APPLICATION
# ============================================================

if __name__ == "__main__":
    main()

