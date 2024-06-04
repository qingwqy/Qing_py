<font face="Helvetica Neue">



### DMRP[^white]

DMRP is a method to calculate the ground state of a quantum many-body system. The idea is to divide the system into two parts, and then calculate the reduced density matrix of one part. The reduced density matrix is then used to calculate the ground state of the whole system. The process is repeated until the whole system is calculated.

Supposed we have a system with $N$ sites in a pure state $\psi$, and we divide the system into two parts, $A$ and $B$. Let $\ket{i}, i = 1,2,\dots,l$ to be a complete set of $A$ and $\ket{j}, i = 1,2,\dots,J$ to be a complete set of $B$. So we can rewrite the state $\psi$ as $ \psi = \sum_{i,j} c_{ij} \ket{i} \ket{j}$.  The reduced density matrix of part $A$ is defined as
$$
\rho_A = \text{Tr}_B \rho = \sum_{j} c_{ij} c_{ij'}^* \ket{i} \bra{j'}.
$$
where $\rho$ is the density matrix of the whole system. 



### Reference 

[^white]:White, Steven R. "Density matrix formulation for quantum renormalization groups." Physical review letters 69.19 (1992): 2863.

</font>