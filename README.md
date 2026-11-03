# PMC Literature Finder

Streamlit aplikacija za pretragu PubMed Central (PMC) baze i izvoz bibliografskih referenci.

## Funkcionalnosti

- Pretraga PMC-a preko NCBI E-utilities (`Bio.Entrez`).
- Izbor broja rezultata (1–100).
- Prikaz PMCID-a, PMID-a (kada je dostupan), naslova, autora, godine, časopisa i DOI-ja.
- Izvoz metapodataka u CSV.
- Generisanje citata u APA 7, Vancouver i Harvard stilu.
- Paket za preuzimanje u ZIP formatu: CSV, RIS, BibTeX, TXT citati u sva tri stila i informacije o upitu.

Aplikacija pretražuje **PMC**, ne celokupni PubMed. Za bibliografske metapodatke koristi PMC ESummary, tako da ne preuzima čitav XML tekst svakog članka niti PDF.

## Instalacija pomoću Conda-e

U direktorijumu projekta:

```bash
conda env create -f environment.yml
conda activate pmc-literature-app
streamlit run app.py
```

## Instalacija u postojećem `bioinfo2026` okruženju

Alternativno, ako želiš da koristiš postojeće okruženje:

```bash
conda activate bioinfo2026
pip install -r requirements.txt
streamlit run app.py
```

## Upotreba

1. Uneti termin ili PMC/Entrez upit, npr. `glioblastoma AND IDH1`.
2. Uneti e-mail adresu koja se šalje NCBI E-utilities servisima.
3. Po želji uneti NCBI API key.
4. Izabrati broj rezultata i stil citiranja.
5. Kliknuti **Pretraži PMC**.
6. Preuzeti CSV, citate ili kompletan ZIP paket.

API key nije obavezan. Nemoj unositi lozinku NCBI naloga. API ključ se koristi samo za zahteve tokom rada aplikacije i nije upisan u projektne fajlove.

## ZIP paket

ZIP sadrži:

- `references.csv` — tabelu bibliografskih metapodataka;
- `references.ris` — RIS za reference menadžere;
- `references.bib` — BibTeX;
- `citations_apa_7.txt`, `citations_vancouver.txt`, `citations_harvard.txt`;
- `search_info.txt` — izvor, upit i broj rezultata.

Automatski formatirani citati su praktična polazna tačka. Pre konačne upotrebe proveriti svaki zapis prema pravilima izabranog časopisa, jer metapodaci i zahtevi stilova mogu varirati.

## Napomena o PMC i pravima korišćenja

Ova verzija preuzima bibliografske metapodatke, ne kopije PDF članaka ili punog teksta. PMC sadrži radove sa različitim licencama; eventualno buduće preuzimanje i redistribuciju punog teksta treba posebno vezati za licencu konkretnog rada.
