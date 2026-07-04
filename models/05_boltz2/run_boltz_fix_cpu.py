import sys
import torch

original_svd = torch.linalg.svd
def safe_svd(*args, **kwargs):
    kwargs.pop("driver", None)  # remove driver arg, not supported with magma
    return original_svd(*args, **kwargs)
torch.linalg.svd = safe_svd

from boltz.main import cli
sys.exit(cli())