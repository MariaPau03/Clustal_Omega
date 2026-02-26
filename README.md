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
```bash
brew install clustal-omega
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
python app.py
# Visit http://localhost:5000
```

**On a server with a custom clustalo path:**
```bash
CLUSTALO_PATH=/opt/bin/clustalo python app.py
```

**Production deployment (gunicorn):**
```bash
pip install gunicorn
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
