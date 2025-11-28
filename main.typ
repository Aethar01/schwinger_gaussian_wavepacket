#import "@preview/physica:0.9.7": *
#let REF = text(fill: red)[REF]
#let appendix(body) = {
  set heading(numbering: "A.1", supplement: "Appendix")
  body
}
#show heading.where(level: 1): set align(center)
#set math.equation(numbering: "(1)")
= Preparation of Gaussian wave packets for the Schwinger model in 1 + 1 dimensions

For the Schwinger model with the $theta$-term in
1 + 1 dimensional U(1) gauge theory, the Lagrangian density is
$
  cal(L)_0 =
  underbrace(
    -1/4 tensor(F, -mu, -nu) tensor(F, mu, nu),
    "Kinetic term\n for the electromagnetic\n gauge field"
  )
  +
  underbrace(
    (g theta) / (4 pi) tensor(epsilon, -mu, -nu) tensor(F, mu, nu),
    "Topological term"
  )
  +
  underbrace(
    i overline(psi) tensor(gamma, mu) (tensor(partial, -mu) + i g tensor(A, mu))psi,
    "Fermion kinetic\n and gauge interaction term"
  )
  underbrace(
    - m overline(psi) psi,
    "Fermion mass term"
  )
$ <l-0>
where
$gamma^0 = sigma^3 = mat(1, 0; 0, -1)$,
$gamma^1 = i sigma^2 = mat(0, -i; i, 0)$ and
$tensor(F, -mu, -nu) = tensor(partial, -mu) tensor(A, -nu) - tensor(partial, -nu) tensor(A, -mu)$.

We can then transform the Lagrangian density by chiral rotation to remove $theta$ from the Gauge Sector:
$
  cal(L) = -1/4 tensor(F, -mu, -nu) tensor(F, mu, nu)
  + i overline(psi) tensor(gamma, mu) (tensor(partial, -mu) + i g tensor(A, -mu))psi

  underbrace(
  - m overline(psi) e^(i theta gamma^5) psi,
  "Modified effective\n mass term"
  )
$
where
$gamma^5 = gamma^0 gamma^1$.

Using an adiabatic method @chakraborty2022 or a VQE (@schwinger-vqe) we can prepare a $ket("vac")$ state for the Schwinger model in 1 + 1 dimensions.

Once we have this state, we need to define some creation operators $c^dagger_n "and" c_n$ which must obey @Chai2025:
$
  {c_n, c^dagger_m} = delta_(n m)
$

@Davoudi2024

#show: appendix
= Additional Information

== Derivation of the Hamiltonian from the Lagrangian density
<schwinger-hamiltonian>
From the Lagrangian density (@l-0) we can derive a Hamiltonian
given in @chakraborty2022. In order to do so we will use the
Legendre transformation $H = integral d^D x cal(H)$ with the
Hamiltonian density in the temporal gauge ($A_0 = 0$):
$
  cal(H) = Pi^A accent(A, dot) + Pi^psi accent(psi, dot) - cal(L)_0
$
where
$Pi^A = (partial cal(L)_0) / (partial accent(A_1, dot)) = tensor(F, 0, 1)$ (in 1 + 1D $tensor(F, 0, 1) = E$),
and
$Pi^psi = (partial cal(L)_0) / (partial accent(psi, dot)) = i overline(psi) tensor(gamma, 0)$

We can now descretize from continuum space to a set of discrete sites $n$ separated by a
lattice spacing $a$: $x -> n a$.
- In temporal gauge, only remaining derivatives are $partial_1 psi(x) -> (psi_(n+1)-psi_n) / a$
- Integrals become sums over lattice sites: $integral d x cal(H)(x) -> a sum_n cal(H)_n$
- Gauge field $A_1$ becomes a link variable between sites $n$ and $n+1$: $A_1 -> U_n = e^(i a g A_1 (x_n))$
- The electric field $E$ becomes an operator $L_n$ between sites $n$ and $n+1$: $E -> L_n = -i (partial) / (partial phi_n)$ so $U_n = e^(-i phi_n)$

  as such, the continuous energy $integral 1/2 E^2 d x -> a sum_n 1/2 E^2_n -> a sum_n 1/2 (-g L_n)^2 = sum_n (g^2 a)/(2) L_n^2$

- Using the staggered fermion formulation, the two-component Dirac spinor $psi (x) -> chi_n$

$
  H = -i sum_n^(N-1) (1/(2 a) - (-1)^n m/2 sin(theta))[chi^dagger_n e^(i
    phi_n) chi_(n+1) - "h.c."] \
  + m cos(theta) sum_n^N (-1)^n chi^dagger_n chi_n
  + (g^2 a)/2 sum_n^(N-1) L^2_n
$

== VQE for the Schwinger model
<schwinger-vqe>
Using the Jordan-Wigner transformations, we can rewrite the fermions in terms of spin variables like:
$
  chi_n = (product_(l<n) - i Z_l) (X_n - i Y_n)/(2)
$

and as in @chakraborty2022:
$
  L_n = L_0 + 1/2 sum_l^N (Z_l + (-1)^l)
$
where we can take $L_0 = 0$ since the model with $(theta, L_0) equiv (theta + 2 pi L_0, 0)$ and redefine $chi_n$ to eliminate $phi_n$: $chi_n -> product_(l<n)[e^(i phi_l)] chi_n$.

Thus, the Hamiltonian becomes:
$
  H = (g^2 a)/4 sum_(n=2)^(N-1) sum_(1 <= k < l <= n) Z_k Z_l \
  + 1/2 sum_(n=1)^(N-1) (1/(2 a) - (-1)^n m/2 sin(theta))[X_n X_(n+1) + Y_n Y_(n+1)] \
  + (m cos(theta)) / 2 sum_(n=1)^N (-1)^n Z_n - (g^2 a)/4 sum_(n=1)^(N-1) (n mod 2) sum_(l=1)^n Z_l
$

Using this Hamiltonian, we can use the VQE to find it's ground state.
1. State preparation
  - Use some ansatz ($U(theta)$) and apply it to some initial state: $U(theta) ket(0)^(times.o N) = ket(psi(theta))$, where $theta$ is a vector of parameters.
2. Expectation value measurement
  - Calculate the expectation value. $E(theta) = braket(psi(theta), H, psi(theta))$
3. Classical optimization
  - The expectation value is then minimized using a classical optimizer. This is done since it is known that the measured energy must be greater than or equal to the true ground state energy.
4. Iterate
  - Repeat steps 2 and 3 until convergence at a minimum.

Our ground state is then defined as $ket(psi(theta_"min"))$.
  

#bibliography(style: "american-physics-society", "bibliography.bib")
