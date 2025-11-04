import torch, torch.nn as nn

class LinearComm(nn.Module):
    def __init__(self, d, gain=0.2):
        super().__init__()
        self.B = nn.Linear(d, d, bias=False)
        with torch.no_grad():
            self.B.weight.copy_(torch.eye(d) * gain)

    def forward(self, Xi, Xj):
        return self.B(Xj - Xi)  # simple, stable f(Δx)

class scQDiffComm(nn.Module):
    def __init__(self, core_model, d):
        super().__init__()
        self.core = core_model
        self.comm = LinearComm(d)

    @torch.no_grad()
    def set_comm_gain(self, g):
        self.comm.B.weight.data.mul_(0.0).add_(torch.eye(self.comm.B.weight.shape[0])*g)

    def drift(self, X, t, Wt):
        u_intra = self.core.drift(X, t)           # reuse your current drift
        Xi, Xj = X.unsqueeze(1), X.unsqueeze(0)   # (N,1,d), (1,N,d)
        Mij = self.comm(Xi.expand_as(Xj), Xj)     # (N,N,d)
        u_comm = torch.einsum('ij,ijd->id', Wt, Mij)
        return u_intra + u_comm
