python -u - << 'PY'
import os
import traceback
from pathlib import Path
from chai_lab.chai1 import run_inference

# Fix memory fragmentation
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

fasta_dir = Path("/home/nehajay/project_new/chai_inputs/")
msa_dir = Path("/home/nehajay/project_new/chai_msas/")
output_root = Path("/home/nehajay/project_new/chai_outputs_28_5/")

print("FASTA:", list(fasta_dir.glob("*.fasta")))

for fasta_file in fasta_dir.glob("*.fasta"):
    print("Running:", fasta_file)
    try:
        run_inference(
            fasta_file=fasta_file,
            output_dir=output_root / fasta_file.stem,
            num_trunk_recycles=1,       # reduced from 3 → saves ~2x memory
            num_diffn_timesteps=100,    # reduced from 200 → faster + less memory
            num_diffn_samples=1,        # reduced from default 5 → big memory saving
            seed=42,
            device="cuda:0",
            use_esm_embeddings=True,
            msa_directory=msa_dir,
        )
        print("DONE:", fasta_file)
    except Exception as e:
        print(f"FAILED: {fasta_file} — {e}")
        traceback.print_exc()
PY