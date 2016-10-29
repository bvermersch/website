We are interested in the ground state properties of the long-range Ising model.
\begin{equation}
H = \sum_{i<j} \frac{J}{(j-i)^3} \sigma_i^x \sigma_j^x +\sum_i \Delta \sigma_i^z
\end{equation}


For $\Delta\to0$, the gap between ground and excited state tend to zero.
Using the standard DMRG approach (without Quantum numbers), the ground state search thus becomes  more and more challenging as  $\Delta\to0$ where the algorithm has to resolve very small energy gaps.
In particular, it becomes extremely hard to access accurately the ground state entropy as a small admixture of the 1st excited state leads to significant errors in the entropy estimation.

This problem can be cured based on the fact that the Hamiltonian conserves the parity of the number of excitations and that ground and first excited states have opposite parity 
(For two particles, these state correspond respectively to $\ket{\uparrow}\ket{\uparrow}-\ket{\downarrow}\ket{\downarrow}$, $\ket{\uparrow}\ket{\downarrow}-\ket{\downarrow}\ket{\uparrow}$).
Using a MPS ansatz with Quantum numbers associated to the parity conservation, the DMRG search will be thus in this case unaffected by the gap closing at $\Delta\to 0$.

\section{Implementation of the site set}
To implement parity conservation with ITensor, we create a new spin set spinhalfparity.h, which is based on spinhalf.h with the crucial differences that the indices describe the number of excitations modulo 2, i.e the parity. This file has to be placed in the folder itensor/mps/sites.

\begin{lstlisting}
inline void SpinHalf::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("S=1/2 ",j),
            Index(nameint("Up ",j),1,Site),QN({1,2}),
            Index(nameint("Dn ",j),1,Site),QN({0,2}));
        }
    }
\end{lstlisting}
and one can use the Quantum number version of the SxSx operator, given that it conserves the parity. This is done by commenting the line  "Op = mixedIQTensor(s,sP);" of the original spinhalf.h class.

\begin{lstlisting}
    if(opname == "Sx")
        {
        //mixedIQTensor call needed here
        //because as an IQTensor, Op would
        //not have a well defined QN flux
        //Op = mixedIQTensor(s,sP);
        Op.set(Up,DnP,+0.5);
        Op.set(Dn,UpP,+0.5);
        }
\end{lstlisting}

\section{Illustration of the performance of the QN approach}

In our code, we now load the new site set

\begin{lstlisting}
#include "itensor/mps/sites/spinhalfparity.h"
#include "itensor/all.h"
\end{lstlisting}

and write our Hamilonian MPO using the autompo module of ITensor:
\begin{lstlisting}
        //Long range ising model
        auto ampo = AutoMPO(sites);
        for(int j = 1; j <= N; ++j)
        {
                ampo += 2.*Delta,"Sz",j;
                for(int k = 1; k <= ncut; ++k)
                        if (j+k<=N)
                                ampo += 4*J*std::pow(k,-alpha),"Sx",j,"Sx",j+k;
        }
        // MPO CONSTRUCTION WITHOUT QN
        auto H = MPO(ampo);
        // MPO WITH QN
        auto HQ = IQMPO(ampo);
\end{lstlisting}

We finally use the dmrg algorithm to compare the efficiency of the parity non-conserving and parity conserving implementations of the MPOs.

\begin{lstlisting}
        //Eigenstate search for non-conserving QN MPO
        Real penalty = 20;
        auto psi = MPS(sites);
        auto en = dmrg(psi,H,sweeps,{"Quiet=",true});
        auto w = std::vector<MPS>(1);
        w.at(0) = psi;
        auto psi1 = MPS(sites);
        auto en1 = dmrg(psi1,H,w,sweeps,{"Quiet=",true,"Weight=",penalty});

        //Eigenstate search conserving QN MPO
        InitState state(sites,"Up");  // Search in the even parity setor
        auto phi = IQMPS(state);
        auto enP = dmrg(phi,HQ,sweeps,{"Quiet=",true});
        state.set(1,"Dn");  // Search in the odd parity sector
        auto phi1 = IQMPS(state);
        auto enP1 = dmrg(phi1,HQ,sweeps,{"Quiet=",true});
     \end{lstlisting}
     
 The gap and the fidelities as a function of the detuning $\Delta$  and for $N=40$ sites are represented in the Figure below showing the advantage of the Quantum number approach.
 
 \includegraphics[width=0.4\columnwidth]{gap}
  \includegraphics[width=0.4\columnwidth]{entropy}
 
Download the new site set.

Download the code.

Download the python scripts to create and plot the data set showed in the Figure.



<span class='article_title>IQTensors: Blocking ITensors by Quantum Number</span>

<span class='article_sig'>Thomas E. Baker&mdash;August 18, 2015</span>

Finding conserved quantities in a Hamiltonian allows for an efficient rewriting of `ITensor`s so that one can identify elements by the numbers denoting these symmettries: the quantum numbers.  Any time that conserved quantities can be identified in a Hamiltonian, it is advantageous to use `IQTensor`s to keep matrix sizes small and computation times faster.  Using the `ITensor`s where an `IQTensor` could be used means that the MPO will carry around a lot of elements that are zero.

In this article, we will detail how to draw diagrams with quantum numbers, identify quantum numbers in a quantum system, and how to program an `IQTensor`.

## How does ITensor use quantum numbers to make calcultions efficient?

Using `ITensor`s for all calculations is a perfectly valid way to write a code.  But there is an inherent inefficiency in doing so.  Let's take a look at an @@S^z@@ operator to see why:

$$
S^z=\begin{pmatrix}
\frac12 & 0\\\\
0 & -\frac12\\\\
\end{pmatrix}
$$

Two of the elements are zero.  It is a waste of memory to then store four values.  In fact, in general, there are a lot of zero entries in an `ITensor`.  `IQTensors` are implemented in our library to take advantage of the fact that we only need to store non-zero entries of an `ITensor`.  In this case, we concentrate on group those non-zero entries by quantum number.

The `IQTensor` will look like a vector divided into quantum numbers:

## Some examples of Quantum Numbers

parity, charge, U(1)

Supeconductiviity

## Diagrams with Quantum Numbers

Diagrams that display the quantum number flux are drawn with arrows.  Fluxes denote the quantum numbers flowing into and out of a block.

<p align="center"><img src="docs/tutorials/IQTensors/iqtensor.png" alt="Diagram" style="width: 400px;"/></p>

    IQMPS psi;
    psi.A(1);

In order to ensure that the flux represented by the arrows in the diagram are correct, we will construct operators that have the corrent flux structure.  This will ensure that our wavefunctions are evaluated in the correct symmetry sector.

The direction of the arrows are a convention.  We want to keep our ket vectors covariant, @@\psi\_\mu@@, so our operators will raise the index, @@\Lambda^{\mu\nu}\psi\_\mu=\psi^\nu@@ on contraction.  The bra vectors are then identified as contravariant vectors.

<p align="center"><img src="docs/tutorials/IQTensors/iqbra.png" alt="Diagram" style="width: 400px;"/></p>

    dag(psi.A(1));

In practice, we don't need to worry about which is which; simply this determines our convention for the direction of arrows. Note that `dag` function reversed the arrow on the physical index.

The link indices (horizontal lines) are chosen from the orthogonality center.  Consider a singlet:

<p align="center"><img src="docs/tutorials/IQTensors/iqnumber1.png" alt="Diagram" style="width: 400px;"/></p>

We can take an SVD of the wavefunction.  The orthogonality center (red diamond) is the source of the quantum numbers (arrows).

<p align="center"><img src="docs/tutorials/IQTensors/iqnumber2.png" alt="Diagram" style="width: 400px;"/></p>

$$
=UDV^\dagger
$$

    IQMPS psi;
    IQTensor U,D,V;
    SVD(psi,U,D,V);

One option is to roll the orthogonality center into the right site.  The arrow then points to the left, naturally. 

<p align="center"><img src="docs/tutorials/IQTensors/iqnumber3.png" alt="Diagram" style="width: 400px;"/></p>

$$
=U(DV^\dagger)=A\Lambda
$$

    auto A = U;//left-normalized tensor
    auto Lambda = D*V;//orthogonality center


We could also have rolled the tensor left:

<p align="center"><img src="docs/tutorials/IQTensors/iqnumber4.png" alt="Diagram" style="width: 400px;"/></p>

$$
=(UD)V^\dagger=\Lambda B
$$

    B = V;//right-normalized
    Lambda = U*D;//orthogonality center

In a less trivial example, the 

<p align="center"><img src="docs/tutorials/IQTensors/iqnumber5.png" alt="Diagram" style="width: 400px;"/></p>

## Operator fluxes 

The flux for the entire Hamiltonian should be zero.  If we had something other than this, we would not have a conservation law! This implies that on the diagram above has a flux of zero on the left and a flux of zero on the right since if this tensor is on either end of the network, this must be the case.  The easiest way to make sure our overall network has a flux of zero is make sure each element has a total flux of zero.

So, the line flowing into the tensor on the left must be zero.

<p align="center"><img src="docs/tutorials/IQTensors/iqflux.png" alt="Diagram" style="width: 400px;"/></p>

The vertical line provides a flux of +1.  It acts like a source of particles.

  <div class="example_clicker">Show Answer</div>

        The missing quantum number on the operator is -1.  Thus, the in quantum number is 1 while the outgoing number is -1.

## Coding IQMPOs


## Not using IQTensors

If necessary, not using the IQTensors for the structure of the MPO is possible.  One simply needs to write the MPO in an approrpriate basis and then use the `dmrg` function.  This invites an increased cost since there are savings, mentioned above, from using quantum numbers.  We strongly recommend using the `IQTensor` feature.




