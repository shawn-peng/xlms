import torch
import torchquad


# Enable GPU support if available and set the floating point precision
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device) # This will print "cuda" if everything work correctly
torchquad.set_up_backend("torch", data_type="float32")
