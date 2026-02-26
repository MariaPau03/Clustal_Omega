# ClustalOmega Web Application

A Flask-based web interface for Multiple Sequence Alignment using **Clustal-Omega**.

## Features

- **4 input modes** (auto-detected):
  - FASTA sequences pasted directly
  - UniProt accession IDs (fetched from `https://www.uniprot.org/uniprot/{id}.fasta`)
  - PDB entry IDs (fetched from `https://www.rcsb.org/fasta/entry/{id}`)
  - FASTA file upload (drag & drop supported)
- **Output formats**: Clustal, FASTA, MSF, PHYLIP, SELEX, Stockholm, Vienna
- **Alignment options**: output format, guide tree iterations, extra clustalo flags
- **Input validation**: format detection, sequence checks, meaningful error messages
- **Downloadable results**

---

## Requirements

- Python 3.8+
- Flask
- Clustal-Omega executable (`clustalo`)

---

## Installation

### 1. Install Clustal-Omega

**Linux (apt):**
```bash
sudo apt install clustalo
```

**Linux (manual):**
```bash
wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64
chmod +x clustalo-1.2.4-Ubuntu-x86_64
sudo mv clustalo-1.2.4-Ubuntu-x86_64 /usr/local/bin/clustalo
```

**macOS (Homebrew):**
Previously create an environment!

```bash
conda create -n clustal_env python=3.10
conda activate clustal_env
```

```bash
conda install -c bioconda clustalo
```

**Verify:**
```bash
clustalo --version
# Expected: Clustal Omega - 1.2.4 (AndreaGiacomo)
```

### 2. Test ClustalOmega from command line

```bash
# Create test input
cat > test.fasta << 'EOF'
>seq1
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKV
>seq2
MHHHHHHGSACPGRHFCCGRGACAVCRWCPGPVQCLDTLTVCPTDCPKERYFCGLYMAMQQLRYNLAF
>seq3
MSLNKIWLQACPGRHFCCGRGACAVCRRCPGPVQCLDTLTVCPTDCPKERYFCGLYMAMQQLRYN
EOF

# Run alignment
clustalo -i test.fasta -o test.aln --outfmt=clustal
cat test.aln
```

### 3. Install Python dependencies

```bash
cd clustal_app/
pip install -r requirements.txt
# Or on systems where pip uses system packages:
pip install --break-system-packages -r requirements.txt
```

### 4. Run the application

```bash
python3 app.py
# Visit http://localhost:5000
```

**On a server with a custom clustalo path:**
```bash
CLUSTALO_PATH=/opt/bin/clustalo python app.py
```

**Production deployment (gunicorn):**
```bash
pip3 install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 app:app
```

---

## Deployment Notes

- Set `CLUSTALO_PATH` env variable if clustalo is not in system PATH
- The `uploads/` and `results/` directories are created automatically
- Result files persist until manually cleaned; consider adding a cron job:
  ```bash
  # Delete result files older than 24h
  find /path/to/clustal_app/results/ -mtime +1 -delete
  ```
- For Apache/Nginx with mod_wsgi, create `wsgi.py`:
  ```python
  from app import app as application
  ```

---

## Project Structure

```
clustal_app/
├── app.py              # Flask application
├── requirements.txt    # Python dependencies
├── templates/
│   └── index.html      # Frontend UI
├── uploads/            # Temporary input files (auto-created)
└── results/            # Output alignment files (auto-created)
```

---

## How to Use the Application

Once the app is running at `http://localhost:5000`, you can submit sequences using any of four input methods. The application **automatically detects** which type you are using — you do not need to select it manually.

---

### Option 1 — Paste FASTA Sequences Directly

If you already have sequences in FASTA format, paste them straight into the text box.

**Format rules:**
- Each sequence must begin with a header line starting with `>`
- The header must be followed immediately by the sequence on the next line(s)
- At least **2 sequences** are required to run an alignment

**Example:**
```
>hemoglobin_alpha_human
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
>hemoglobin_beta_human
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH
```

> **Tip:** If you copy a sequence from a database webpage, make sure the `>header` line is included. A common mistake is copying only the sequence letters without the header — the app will detect this and show a clear error.

---

### Option 2 — UniProt Accession IDs

[UniProt](https://www.uniprot.org) is the primary database for protein sequences and functional annotation. The app will fetch each sequence automatically using the UniProt REST API.

**How to find a UniProt ID:**

1. Go to [https://www.uniprot.org](https://www.uniprot.org)
2. Search for your protein by name, gene, or organism (e.g. `hemoglobin alpha human`)
3. Click on a result — the accession ID is the short code at the top (e.g. **P69905**)
4. It also appears in the URL: `https://www.uniprot.org/uniprot/P69905`

**Format:** Enter one ID per line (or space/comma separated). IDs are 6–10 alphanumeric characters.

```
P69905
P68871
P02042
P01942
```

**What the app does:** fetches `https://www.uniprot.org/uniprot/{ID}.fasta` for each ID, combines the sequences, and runs the alignment. If one ID is invalid, the others still proceed and you receive a warning.

> **Tip:** UniProt has two sections — **Swiss-Prot** (manually reviewed, higher quality, e.g. `P69905`) and **TrEMBL** (automatically annotated, e.g. `A0A000`). Prefer Swiss-Prot entries when available — look for the gold star ⭐ badge on the search results page.

> **Note:** UniProt and PDB always return **protein** sequences. If you have DNA or RNA selected in the Sequence Type option, the app will automatically override it to Protein and show a warning.

---

### Option 3 — PDB Entry IDs

The [RCSB Protein Data Bank](https://www.rcsb.org) (PDB) stores 3D structures of proteins, nucleic acids, and complexes. Each entry has a sequence that the app can fetch for alignment.

**How to find a PDB ID:**

1. Go to [https://www.rcsb.org](https://www.rcsb.org)
2. Search for your protein (e.g. `subtilisin inhibitor` or `1SBI`)
3. The PDB ID is the **4-character code** shown at the top of the entry page (e.g. **1SBI**)
4. It also appears in the URL: `https://www.rcsb.org/structure/1SBI`

**Format:** Enter one 4-character ID per line.

```
1HHO
2HHB
1A3N
1GZX
```

**What the app does:** fetches `https://www.rcsb.org/fasta/entry/{ID}` for each entry. Note that PDB entries often contain **multiple chains** (e.g. chain A, chain B) — all chains from the entry are included in the FASTA returned by RCSB. If you only want a specific chain, use the FASTA sequence directly (Option 1) after downloading it from the PDB page.

> **Tip:** You can also enter PDB IDs with a chain letter appended (e.g. `1SBIA`) — the app will automatically strip the chain letter and fetch the full entry.

> **Tip:** To compare the same protein from different organisms or crystal forms, collect their PDB IDs from the RCSB search results and paste them all at once.

---

### Option 4 — Upload a FASTA File

Use this option when you have a `.fasta` file on your computer — for example, downloaded from a database or generated by another tool.

**How to use:**
1. Click the **File Upload** tab at the top of the input panel
2. Drag and drop your file onto the upload zone, or click to browse
3. Accepted formats: `.fasta` `.fa` `.fas` `.txt` `.seq`
4. Maximum file size: 16 MB

**How to download a FASTA file from UniProt:**
1. Open a UniProt entry (e.g. `https://www.uniprot.org/uniprot/P69905`)
2. Click **Download** → select format **FASTA (canonical)** → click **Download**

**How to download a FASTA file from RCSB PDB:**
1. Open a PDB entry (e.g. `https://www.rcsb.org/structure/1HHO`)
2. Click **Download Files** → select **FASTA Sequence**

**How to download from NCBI:**
1. Go to [https://www.ncbi.nlm.nih.gov/protein](https://www.ncbi.nlm.nih.gov/protein)
2. Search for your protein
3. Click **Send to** → **File** → format **FASTA** → **Create File**

---

### Alignment Options

Once your input is ready, configure the alignment in the **Alignment Options** panel:

| Option | What it does |
|---|---|
| **Sequence Type** | Select Protein, DNA, or RNA. This sets the `--seqtype` flag in clustalo and determines which characters are considered valid. |
| **Output Format** | Choose how the alignment result is formatted. **Clustal** (.aln) is the default and shows a conservation line below the sequences. Other options include FASTA, PHYLIP, Stockholm, and more. |
| **Guide Tree Iterations** | Set to 0 for a standard single-pass alignment. Increasing to 1–5 causes clustalo to refine the alignment by rebuilding the guide tree — useful for distantly related sequences, at the cost of extra runtime. |
| **Extra Options** | Advanced: pass any additional `clustalo` command-line flags directly, e.g. `--threads=4` to use multiple CPU cores. |

---

### Reading the Results

After a successful run:

- The **stats bar** shows how many sequences were aligned, their length range, and which input type and format was used
- The **Formatted tab** displays the alignment with colour-coded conservation symbols (Clustal format only):
  - `*` — all sequences have the same residue at this position
  - `:` — residues at this position are strongly similar
  - `.` — residues at this position are weakly similar
  - ` ` (space) — no conservation
- The **Raw tab** shows the plain output text, which you can copy directly
- The **Download button** saves the alignment file to your computer with the correct extension for the chosen format

---

## Input Format Detection

The app automatically detects input type without user intervention:

| Pattern | Detected As |
|---------|-------------|
| Lines starting with `>` | FASTA sequences |
| 6-10 char alphanumeric tokens matching UniProt regex | UniProt IDs |
| 4-char tokens matching `[0-9][A-Z0-9]{3}` | PDB IDs |

---

## Error Handling

- Empty input → clear message
- Unknown format → descriptive error
- < 2 sequences → minimum requirement message
- Invalid amino acid characters → reports specific characters and sequence
- UniProt/PDB fetch failures → per-ID errors, partial success still proceeds
- ClustalOmega not found → informative message with setup instructions
- ClustalOmega runtime errors → stderr output shown to user
- Unsafe extra options (shell injection chars) → rejected with message
