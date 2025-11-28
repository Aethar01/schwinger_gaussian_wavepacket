#import "@preview/physica:0.9.7": *
#let REF = text(fill: red)[REF]
#let appendix(body) = {
  set heading(numbering: "A.1", supplement: "Appendix")
  body
}
#show heading.where(level: 1): set align(center)
#set math.equation(numbering: "(1)")
#let nonum(eq) = math.equation(block: true, numbering: none, eq)
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

Once we have this state, we need to define some creation operators $chi^dagger_n "and" chi_n$ which must obey @Chai2025:
$
  {chi_n, chi^dagger_m} = delta_(n m)
$ <anticommutator>

From the Jordan-Wigner transformation, these creation operators could be represented as:
$
  chi_n = (product_(l<n) - i Z_l) (X_n - i Y_n)/2 \
  chi^dagger_n = (product_(l<n) i Z_l) (X_n + i Y_n)/2
$

Which obeys @anticommutator and is shown in @anticommutator-proof.

Now that we have the creation and annihilation operators, we can create a pair of operators that both act on $ket("vac")$ to create a Gaussian wave packet of either fermions or anti-fermions.

$
  C^dagger = sum_(n=1)^N f^c_n chi^dagger_n, space D^dagger = sum_(n=1)^N f^d_n chi_n
$

where $f^(c(d))_n$ is the Gaussian coefficient at site $n$.

$
  f^c_n (n_0, sigma, k_0) = A^c
  dot 
  exp[-(n - n_0)^2 / (4 sigma^2)] dot exp[i k_0 a(n - n_0)]
  dot
  sqrt((m_"eff" + w_(k_0))/w_(k_0))
  dot
  (Pi_(n 0) + v_(k_0) Pi_(n 1))

  \
  f^d_n (n_0, sigma, k_0) = A^d
  dot
  exp[-(n - n_0)^2 / (4 sigma^2)]
  dot
  exp[-i k_0 a(n - n_0)]
  dot
  sqrt((m_"eff" + w_(k_0))/w_(k_0))
  dot
  (Pi_(n 1) + v_(k_0) Pi_(n 0))
$

Where A is the normalization constant for the condition $sum_(n=1)^N abs(f_n)^2 = 1$, $n_0$ is the center position of the wave packet, $sigma$ is the standard deviation, $k_0$ is the average momentum of the wave packet, $a$ is the lattice spacing, $m_"eff" = m cos(theta)$ is the effective mass, $w_(k_0) = sqrt((m cos theta)^2 + sin^2 k_0)$ and $v_(k_0) = (sin k_0) / (m cos theta + w_(k_0))$ are the group and phase velocity and $Pi_(n l)$ are the projection operators defined as in @Chai2025:
$
  Pi_(n l) = (1 + (-1)^(n+l))/2, space l in {0, 1}
$

Deriving $A^c$:
$
  sum_(n=1)^N abs(
    A^c
    dot
    exp[-(n - n_0)^2 / (4 sigma^2)]
    dot
    sqrt((m_"eff" + w_(k_0))/w_(k_0))
    dot
    (Pi_(n 0) + v_(k_0) Pi_(n 1))
    dot
    underbrace(
      exp[i k_0 a(n - n_0)],
      "magnitude of 1"
    )
  )^2 = 1\

  (A^c)^2
  dot
  sum_(n=1)^N exp[-(n - n_0)^2 / (2 sigma^2)]
  dot
  (m_"eff" + w_(k_0))/w_(k_0)
  dot
  (Pi_(n 0) + v_(k_0) Pi_(n 1))^2
  = 1\
  
  A^c (N, sigma, n_0) = (
    sum_(n=1)^N exp[-(n - n_0)^2 / (2 sigma^2)]
    dot
    (m_"eff" + w_(k_0))/w_(k_0)
    dot
    (Pi_(n 0) + v_(k_0) Pi_(n 1))^2
  )^(-1/2)
$
Which then follows that $A^d$ is:
$
  A^d (N, sigma, n_0) = (
    sum_(n=1)^N exp[-(n - n_0)^2 / (2 sigma^2)]
    dot
    (m_"eff" + w_(k_0))/w_(k_0)
    dot
    (Pi_(n 1) + v_(k_0) Pi_(n 0))^2
  )^(-1/2)
$

We can now use these operators to create a fermion, anti-fermion wave packet on top of the vacuum state:
$
  ket(psi(0)) = D^dagger C^dagger ket("vac")
$

// @Davoudi2024

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
#nonum[
  $
    chi_n = (product_(l<n) - i Z_l) (X_n - i Y_n)/(2)
  $
]

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

== Anti-commutation relation for the creation and annihilation operators
<anticommutator-proof>

Our creation and annihilation operators are:
#nonum[
  $
    chi_n = (product_(l<n) - i Z_l) (X_n - i Y_n)/2 \
    chi^dagger_n = (product_(l<n) i Z_l) (X_n + i Y_n)/2
  $
]

And we require that:
#nonum[
  $
    {chi_n, chi^dagger_m} = delta_(n m)
  $
]

Let's simplify the expression slightly, defining $P_n$ and $P^dagger_n$:
$
  P_n &= product_(l<n) - i Z_l\
  P^dagger_m &= product_(l<m) i Z_l
$

Substituting in:
$
  {chi_n, chi^dagger_m} = chi_n chi^dagger_m + chi^dagger_m chi_n\
  = (P_n (X_n - i Y_n)/2) (P^dagger_m (X_m + i Y_m)/2) + (P^dagger_m (X_m + i Y_m)/2) (P_n (X_n - i Y_n)/2)\

  "and since:"\

  {Z_n, (X_n - i Y_n)/2} = mat(1,0;0,-1) mat(0,0;1,0) + mat(0,0;1,0) mat(1,0;0,-1) = mat(0,0;-1,0) + mat(0,0;1,0) = 0\

  "we can rewrite as:"\

  {chi_n, chi^dagger_m} = 1/4 [(P_n P^dagger_m) (X_n - i Y_n) (X_m + i Y_m) + (P^dagger_m P_n) (X_m + i Y_m) (X_n - i Y_n)]\

  "and since:"\

  i Z (i Z)^dagger = (i Z)^dagger i Z = - (i Z) (i Z) = -i^2 Z^2 = I\

  "for" n=m "we can rewrite as:"\
  
  {chi_n, chi^dagger_n} = 1/4 [(X_n - i Y_n) (X_n + i Y_n) + (X_n + i Y_n) (X_n - i Y_n)]\
  = 1/4 [X_n^2 + Y_n^2 + i[X_n, Y_n] + X_n^2 + Y_n^2 - i[X_n, Y_n]]\
  = 1/4 [2 X_n^2 + 2 Y_n^2]\
  = I_(n n) = 1\
$ <nem>
$
  "and for" n != m "assuming" n < m "we can rewrite as:"\

  {chi_n, chi^dagger_m} = (P_n (X_n - i Y_n)/2) (P^dagger_m (X_m + i Y_m)/2) + (P^dagger_m (X_m + i Y_m)/2) (P_n (X_n - i Y_n)/2)\

  "we can the split P_m into two parts:"\

  P^dagger_m = P^dagger_n (i Z_n) P'^dagger_(n,m)\

  "where" P^dagger_n "acts on sites" l<n "and" P'^dagger_(n,m) "acts on sites" n<l<m\
  "if" F_n = (X_n - i Y_n)/2, F^dagger_m = (X_m + i Y_m)/2 "then:"\
  
  {chi_n, chi^dagger_m} = 1/4 [(P_n P^dagger_n) F_n (i Z_n) P'^dagger_(n,m) F^dagger_m + (P^dagger_n P_n) (i Z_n) P'^dagger_(n,m) F^dagger_m F_n]\
  
  "since" P_n P^dagger_n = P^dagger_n P_n = I ":"\

  {chi_n, chi^dagger_m} = i/4 [F_n Z_n (P'^dagger_(n,m) F^dagger_m) + Z_n F_n (P'^dagger_(n,m) F^dagger_m)]\

  "and as in" #[@nem], {Z_n, F_n} = 0\

  {chi_n, chi^dagger_m} = i/4 [(-Z_n F_n) (P'^dagger_(n,m) F^dagger_m) + Z_n F_n (P'^dagger_(n,m) F^dagger_m)]\
  = 0
$

Therefore, ${chi_n, chi^dagger_m} = delta_(n m)$.

  
// 1. For $n = m$ the $Z$ strings cancel and the resulting expression is:
//   $
//     {chi_n, chi^dagger_n} = 1/2 {(X_n - i Y_n), (X_n + i Y_n)} \
//   $
// 2. We can then solve using Pauli matrices:
//   $
//     (X - i Y)/2 = mat(0,1;0,0),
//     (X + i Y)/2 = mat(0,0;0,1)\
//     1/2 ((X - i Y) (X + i Y) + (X + i Y) (X - i Y)) = mat(1,0;0,1)\
//     "therefore" chi_n chi^dagger_n + chi^dagger_n chi_n = 1\
//   $

#bibliography(style: "american-physics-society", "bibliography.bib")
