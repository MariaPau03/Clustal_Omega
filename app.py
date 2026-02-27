#!/usr/bin/env python3
"""
ClustalOmega Web Application
Flask-based interface for multiple sequence alignment
"""

import os
import re
import uuid
import subprocess
import tempfile
import requests
import time
from flask import Flask, render_template, request, jsonify, send_file, abort
from werkzeug.utils import secure_filename



app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max upload
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'uploads')
app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(__file__), 'results')

# Ensure directories exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

# Path to clustal-omega executable (adjust for server)
CLUSTALO_PATH = os.environ.get('CLUSTALO_PATH', 'clustalo')

# Allowed file extensions
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fas', 'txt', 'seq'}

# Output format options
OUTPUT_FORMATS = {
    'clustal': 'Clustal (.aln)',
    'fasta': 'FASTA (.fasta)',
    'msf': 'MSF (.msf)',
    'phylip': 'PHYLIP (.phy)',
    'selex': 'SELEX (.slx)',
    'stockholm': 'Stockholm (.sto)',
    'vienna': 'Vienna (.vienna)',
}

FORMAT_EXTENSIONS = {
    'clustal': 'aln',
    'fasta': 'fasta',
    'msf': 'msf',
    'phylip': 'phy',
    'selex': 'slx',
    'stockholm': 'sto',
    'vienna': 'vienna',
}

# Sequence type options (value → clustalo --seqtype argument)
SEQUENCE_TYPES = {
    'protein': 'Protein',
    'dna':     'DNA',
    'rna':     'RNA',
}

# Valid residue characters per sequence type (IUPAC, case-insensitive)
SEQ_TYPE_CHARS = {
    'protein': re.compile(r'[^ACDEFGHIKLMNPQRSTVWYXBZUJacdefghiklmnpqrstvwyxbzuj*\-]'),
    'dna':     re.compile(r'[^ACGTURYSWKMBDHVNacgturyswkmbdhvn\-]'),
    'rna':     re.compile(r'[^ACGURYSWKMBDHVNacguryswkmbdhvn\-]'),
}

SEQ_TYPE_LABELS = {
    'protein': 'Protein',
    'dna':     'DNA',
    'rna':     'RNA',
}

# ─── Sequence Detection & Validation ─────────────────────────────────────────

def detect_input_type(text):
    """Auto-detect whether input is FASTA sequences, UniProt IDs, or PDB IDs."""
    text = text.strip()
    if not text:
        return None, "Input is empty."

    lines = [l.strip() for l in text.splitlines() if l.strip()]

    # Check if it looks like FASTA
    if lines[0].startswith('>'):
        return 'fasta', None

    # Try to identify IDs (one per line or space-separated)
    tokens = []
    for line in lines:
        tokens.extend(line.split())

    if not tokens:
        return None, "No identifiable input found."

    # UniProt pattern: 6-10 alphanumeric chars (e.g., P12345, Q9Y6K9, A0A000AAA0)
    uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', re.IGNORECASE)
    # PDB pattern: 4-char alphanumeric (e.g., 1ABC, 2XYZ)
    pdb_pattern = re.compile(r'^[0-9][A-Z0-9]{3}$', re.IGNORECASE)

    uniprot_matches = sum(1 for t in tokens if uniprot_pattern.match(t))
    pdb_matches = sum(1 for t in tokens if pdb_pattern.match(t))

    if uniprot_matches > 0 and uniprot_matches >= pdb_matches:
        return 'uniprot', None
    elif pdb_matches > 0:
        return 'pdb', None
    else:
        # Could still be UniProt IDs with unusual format - try broader check
        # UniProt IDs are 6-10 chars
        broad_uniprot = sum(1 for t in tokens if re.match(r'^[A-Z0-9]{6,10}$', t, re.IGNORECASE))
        if broad_uniprot == len(tokens):
            return 'uniprot', None
        return None, f"Unrecognized input format. Could not identify as FASTA, UniProt IDs, or PDB IDs. Got: {tokens[:3]}"


def validate_fasta(fasta_text, seq_type='protein'):
    """Validate FASTA format for the given sequence type and return (sequences_dict, error)."""
    sequences = {}
    current_id = None
    current_seq = []

    bad_char_re = SEQ_TYPE_CHARS.get(seq_type, SEQ_TYPE_CHARS['protein'])
    type_label  = SEQ_TYPE_LABELS.get(seq_type, 'Protein')

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_id:
                seq = ''.join(current_seq)
                if not seq:
                    return None, f"Sequence '{current_id}' has no residues."
                sequences[current_id] = seq
            current_id = line[1:].strip()
            if not current_id:
                return None, "Found a '>' header with no sequence ID."
            current_seq = []
        else:
            if current_id is None:
                return None, "Sequence data found before any FASTA header ('>...')."
            cleaned = re.sub(r'\s', '', line)
            invalid = bad_char_re.sub('', cleaned)   # characters NOT in the valid set
            # Re-find what was removed vs original
            found_invalid = bad_char_re.findall(cleaned)
            if found_invalid:
                bad_sample = ''.join(dict.fromkeys(found_invalid))[:10]
                return None, (
                    f"Invalid {type_label} characters in sequence '{current_id}': "
                    f"'{bad_sample}'. Check that the correct sequence type is selected."
                )
            current_seq.append(cleaned)

    if current_id:
        seq = ''.join(current_seq)
        if not seq:
            return None, f"Sequence '{current_id}' has no residues."
        sequences[current_id] = seq

    if len(sequences) < 2:
        return None, f"At least 2 sequences are required for alignment. Found: {len(sequences)}."

    return sequences, None


def fetch_uniprot(uid):
    """Fetch FASTA sequence from UniProt."""
    uid = uid.strip().upper()
    url = f"https://www.uniprot.org/uniprot/{uid}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200 and resp.text.startswith('>'):
            return resp.text, None
        elif resp.status_code == 404:
            return None, f"UniProt ID '{uid}' not found."
        else:
            return None, f"UniProt returned status {resp.status_code} for '{uid}'."
    except requests.exceptions.Timeout:
        return None, f"Timeout fetching UniProt ID '{uid}'."
    except Exception as e:
        return None, f"Error fetching UniProt ID '{uid}': {str(e)}"


def fetch_pdb(pid):
    """Fetch FASTA sequence from RCSB PDB."""
    pid = pid.strip().upper()
    url = f"https://www.rcsb.org/fasta/entry/{pid}"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200 and resp.text.strip().startswith('>'):
            return resp.text, None
        elif resp.status_code == 404:
            return None, f"PDB ID '{pid}' not found."
        else:
            return None, f"RCSB returned status {resp.status_code} for '{pid}'."
    except requests.exceptions.Timeout:
        return None, f"Timeout fetching PDB ID '{pid}'."
    except Exception as e:
        return None, f"Error fetching PDB ID '{pid}': {str(e)}"


def fetch_sequences_from_ids(ids, id_type):
    """Fetch multiple sequences from UniProt or PDB IDs."""
    combined_fasta = ""
    errors = []
    fetched = 0

    for uid in ids:
        uid = uid.strip()
        if not uid:
            continue
        if id_type == 'uniprot':
            fasta, err = fetch_uniprot(uid)
        else:
            fasta, err = fetch_pdb(uid)

        if err:
            errors.append(err)
        else:
            combined_fasta += fasta.strip() + "\n"
            fetched += 1

    if errors and fetched == 0:
        return None, "Failed to fetch any sequences:\n" + "\n".join(errors)
    if fetched < 2:
        msg = f"Only {fetched} sequence(s) fetched successfully. Need at least 2."
        if errors:
            msg += "\nErrors:\n" + "\n".join(errors)
        return None, msg

    warnings = errors  # partial failures become warnings
    return combined_fasta, warnings if warnings else None


# ─── ClustalOmega Runner ──────────────────────────────────────────────────────

def check_clustalo():
    """Check if clustalo is available."""
    try:
        result = subprocess.run(
            [CLUSTALO_PATH, '--version'],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0, result.stdout.strip()
    except FileNotFoundError:
        return False, None
    except Exception:
        return False, None


def run_clustalo(fasta_text, out_format='clustal', seq_type='protein', extra_opts='', iterations=0):
    """
    Run Clustal-Omega and return (result_text, result_path, error).
    """
    job_id = str(uuid.uuid4())[:8]
    input_path = os.path.join(app.config['UPLOAD_FOLDER'], f"input_{job_id}.fasta")
    ext = FORMAT_EXTENSIONS.get(out_format, 'aln')
    output_path = os.path.join(app.config['RESULTS_FOLDER'], f"result_{job_id}.{ext}")

    # Write input
    with open(input_path, 'w') as f:
        f.write(fasta_text)

    # Map internal key to clustalo --seqtype value
    seqtype_arg = SEQUENCE_TYPES.get(seq_type, 'Protein')

    # Build command
    cmd = [
        CLUSTALO_PATH,
        '-i', input_path,
        '-o', output_path,
        '--outfmt', out_format,
        '--seqtype', seqtype_arg,
        '--force',
    ]

    if iterations > 0:
        cmd += ['--iter', str(iterations)]

    # Parse extra options safely
    if extra_opts.strip():
        # Split carefully; prevent injection
        safe_opts = extra_opts.strip()
        # Only allow safe characters
        if re.search(r'[;&|`$<>]', safe_opts):
            return None, None, "Extra options contain unsafe characters."
        import shlex
        try:
            parsed = shlex.split(safe_opts)
            cmd += parsed
        except ValueError as e:
            return None, None, f"Invalid extra options: {e}"

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )

        # Clean up input file
        try:
            os.remove(input_path)
        except Exception:
            pass

        if result.returncode != 0:
            err_msg = result.stderr.strip() or result.stdout.strip()
            return None, None, f"ClustalOmega error (code {result.returncode}):\n{err_msg}"

        if not os.path.exists(output_path):
            return None, None, "ClustalOmega ran but produced no output file."

        with open(output_path, 'r') as f:
            output_text = f.read()

        return output_text, output_path, None

    except subprocess.TimeoutExpired:
        return None, None, "ClustalOmega timed out after 120 seconds. Your input may be too large."
    except FileNotFoundError:
        return None, None, f"ClustalOmega executable not found at '{CLUSTALO_PATH}'. Please ensure it is installed and in PATH."
    except Exception as e:
        return None, None, f"Unexpected error running ClustalOmega: {str(e)}"


# ─── Routes ───────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    available, version = check_clustalo()
    return render_template('index.html',
                           output_formats=OUTPUT_FORMATS,
                           sequence_types=SEQUENCE_TYPES,
                           clustalo_available=available,
                           clustalo_version=version)


@app.route('/align', methods=['POST'])
def align():
    """Main alignment endpoint."""
    errors = []
    warnings = []

    # ── Gather input ──
    input_mode = request.form.get('input_mode', 'text')
    out_format = request.form.get('out_format', 'clustal')
    seq_type   = request.form.get('seq_type', 'protein').lower()
    extra_opts = request.form.get('extra_opts', '')
    iterations = int(request.form.get('iterations', 0))

    if out_format not in OUTPUT_FORMATS:
        return jsonify({'success': False, 'error': f"Unknown output format: '{out_format}'"}), 400

    if seq_type not in SEQUENCE_TYPES:
        return jsonify({'success': False, 'error': f"Unknown sequence type: '{seq_type}'. Choose protein, dna, or rna."}), 400

    fasta_text = None

    if input_mode == 'file':
        # File upload
        if 'fasta_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file was uploaded.'}), 400
        f = request.files['fasta_file']
        if f.filename == '':
            return jsonify({'success': False, 'error': 'No file selected.'}), 400
        ext = f.filename.rsplit('.', 1)[-1].lower() if '.' in f.filename else ''
        if ext not in ALLOWED_EXTENSIONS:
            return jsonify({'success': False, 'error': f"File type '.{ext}' not allowed. Use FASTA format (.fasta, .fa, .fas, .txt)."}), 400
        try:
            fasta_text = f.read().decode('utf-8', errors='replace')
        except Exception as e:
            return jsonify({'success': False, 'error': f"Could not read uploaded file: {e}"}), 400

        input_type = 'fasta'

    else:
        # Text input
        raw_text = request.form.get('sequences', '').strip()
        if not raw_text:
            return jsonify({'success': False, 'error': 'Input is empty. Please provide sequences or IDs.'}), 400

        input_type, det_error = detect_input_type(raw_text)
        if det_error:
            return jsonify({'success': False, 'error': det_error}), 400

        if input_type in ('uniprot', 'pdb'):
            # UniProt and PDB always return protein sequences
            if seq_type in ('dna', 'rna'):
                warnings.append(
                    f"UniProt and PDB entries contain protein sequences. "
                    f"Sequence type has been overridden from '{seq_type.upper()}' to 'Protein'."
                )
                seq_type = 'protein'

            # Extract IDs
            ids = re.split(r'[\s,;]+', raw_text)
            ids = [i.strip() for i in ids if i.strip()]
            if len(ids) < 2:
                return jsonify({'success': False, 'error': f"Need at least 2 {input_type.upper()} IDs. Got {len(ids)}."}), 400

            fetched, fetch_warn = fetch_sequences_from_ids(ids, input_type)
            if not fetched:
                return jsonify({'success': False, 'error': fetch_warn}), 400
            if fetch_warn:
                warnings.extend(fetch_warn if isinstance(fetch_warn, list) else [fetch_warn])
            fasta_text = fetched

        else:
            fasta_text = raw_text

    # ── Validate FASTA ──
    sequences, val_error = validate_fasta(fasta_text, seq_type)
    if val_error:
        return jsonify({'success': False, 'error': f"Sequence validation error: {val_error}"}), 400

    seq_count = len(sequences)

    # ── Run ClustalOmega ──
    result_text, result_path, run_error = run_clustalo(
        fasta_text, out_format, seq_type, extra_opts, iterations
    )

    if run_error:
        return jsonify({'success': False, 'error': run_error}), 500

    # ── Build stats ──
    lengths = [len(v) for v in sequences.values()]
    stats = {
        'sequences': seq_count,
        'min_length': min(lengths),
        'max_length': max(lengths),
        'avg_length': round(sum(lengths) / len(lengths)),
        'format': OUTPUT_FORMATS.get(out_format, out_format),
        'seq_type': SEQ_TYPE_LABELS.get(seq_type, seq_type.capitalize()),
        'result_file': os.path.basename(result_path),
    }

    return jsonify({
        'success': True,
        'result': result_text,
        'stats': stats,
        'warnings': warnings,
        'input_type': input_type,
        'out_format': out_format,
        'seq_type': seq_type,
        'result_file': os.path.basename(result_path),
    })


@app.route('/download/<filename>')
def download(filename):
    """Download a result file."""
    # Security: only allow safe filenames from our results folder
    if not re.match(r'^result_[a-f0-9]{8}\.\w+$', filename):
        abort(403)
    filepath = os.path.join(app.config['RESULTS_FOLDER'], filename)
    if not os.path.exists(filepath):
        abort(404)
    return send_file(filepath, as_attachment=True, download_name=filename)


@app.route('/status')
def status():
    """Check ClustalOmega availability."""
    available, version = check_clustalo()
    return jsonify({'available': available, 'version': version})


# Ensure these imports are at the top or here
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.exceptions import NotFound

PREFIX = '/u316755/clustal'

# This must be the part uWSGI loads
hostedApp = Flask(__name__)
hostedApp.wsgi_app = DispatcherMiddleware(NotFound(), {
    PREFIX: app  # 'app' is your original Flask(name) variable
})

if __name__ == "__main__":
    app.run()
