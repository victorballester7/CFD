#import "math_commands.typ": *
#import "@preview/physica:0.9.2": *
#show: thmrules

#show list: it => pad(left: 2.5em, {
    set list(marker: style(styles => place(dx: -measure([•], styles).width, [•])))
    it
  })

// Label equations and references ------------
#show math.equation:it => {
  if it.has("label") {
    math.equation(block: true, numbering: "(1)", supplement: [Eq.], it)
  } else {
    it
  }
}

#show ref: it => {
  let el = it.element
  if el != none and el.func() == math.equation {
    link(
      el.location(),
      numbering("Eq. (1)", counter(math.equation).at(el.location()).at(0) + 1),
    )
  } else {
    it
  }
}

// style
#set heading(numbering: "1.1.") // set heading numbering
#set page(paper: "a4", margin: auto,numbering: "1")



#align(center, text(25pt)[
  *Tutorial 4*
])
#align(center, text(15pt)[
  Víctor Ballester Ribó\ \
  Numerical Methods for Fluid Dynamics\ 
  M2 - Applied and Theoretical Mathematics\ 
  Université Paris-Dauphine\ 
  March 2024\ 
])

#exercice[
  Let us consider the heat equation in an open bounded domain $Omega subset RR^n$, over a time interval $[0, T]$, with a homogeneous Dirichlet boundary condition
  $
    cases(
      L(u) = pdv(u, t) - kappa pdv(u, x,2)quad  &"if" (x,t) in Omega times (0, T],
      u = 0 quad &"if" diff Omega times (0, T],
      u(x,0) = u_0(x) quad &"if" x in Omega
    )
  $ <eq:heat>
  We remind the _maximum principle_ for this continuous problem, i.e.
  $
    min lr((0,inf_(x in Omega) u_0(x))) <= u(x,t) <= max lr((0,sup_(x in Omega) u_0(x))), quad forall (x,t) in Omega times (0,T]
  $

  1. We consider @eq:heat, and want to show that the explicit Euler scheme
    $
      L_h (u) = (u_j^(n+1)-u_j^n)/(Delta t) - kappa (u_(j+1)^(n+1) - 2 u_j^(n+1) + u_(j-1)^(n+1))/(Delta x^2) = 0
    $
    satisfies, under some conditions (to be determined) the discrete maximum principle, i.e.
    $
      min lr((0,min_(j in [1,N]) u_j^0)) <= u_j^n <= max lr((0,max_(j in [1,N]) u_j^0)), quad forall j in [1,N] times [0,M]
    $
    for all initial condition $u_0$. In other words, if the initial data is bounded by two constants $m <= 0 <= M$, 
    $
      m <= u_j^0 <= M, quad forall j in [1,N]
    $
    then 
    $
      m <= u_j^n <= M, quad forall j in [1,N] and n >= 0
    $
    #solution[
      We proceed by induction. The case $n=0$, it's already given so we assume that the result holds for $n$ and we prove it for $n+1$. Let $lambda := (Delta t)/(Delta x)^2 kappa$ and suppose we have $lambda < 1/2$. Then:
      $
        u_j^n + lambda (u_(j+1)^n - 2 u_j^n + u_(j-1)^n)&= u_j^(n+1)= u_j^n + lambda (u_(j+1)^n - 2 u_j^n + u_(j-1)^n) \
        (1 - 2 lambda) u_j^n + lambda (u_(j+1)^n + u_(j-1)^n)&= u_j^(n+1)= (1 - 2 lambda) u_j^n + lambda (u_(j+1)^n + u_(j-1)^n)\
        (1 - 2 lambda) m + lambda (m + m)&<= u_j^(n+1)<= (1 - 2 lambda) M + lambda (M + M)\
        m&<= u_j^(n+1)<= M
      $
    ]
  2. Let us define the following discrete norms, for all $1 <= p <= oo$,
    $
      ||u^n||_p = lr((sum_(j=1)^N |u_j^n|^p Delta x))^(1/p)
    $
    and $||u^n||_oo = max_(j in [1,N]) |u_j^n|$. A numerical scheme is said to be stable in the $||dot.c||$ norm if there exists a constant $K>0$, independent of the $Delta x$ and $Delta t$, such that for all initial conditions $u_0$,
    $
      ||u^n|| <= K ||u^0||, quad forall n >= 0
    $ <eq:stability>
    Deduce the stability condition in the $L^oo$ norm for the explicit Euler scheme. Compare to the criterion derived during the lecture for stability in the $L^2$ norm.
    #solution[
      We have that
      $
        ||u^(n+1)||_oo &= max_(j in [1,N]) |u_j^(n+1)| \
        &<= |1 - 2 lambda| max_(j in [1,N]) |u_j^n| + lambda max_(j in [0,N-1]) |u_(j)^n| + lambda max_(j in [2,N+1])|u_(j+1)^n| \
        &= |1 - 2 lambda| ||u^n||_oo + 2 lambda ||u^n||_oo
      $
      Now if $lambda <= 1/2$, we get $||u^(n+1)||_oo <= ||u^n||_oo$, which iteratively implies $||u^n||_oo <= ||u^0||_oo$ (with $K=1$). Now assume $lambda > 1/2$. Assume we have an initial condition $u_0$ such that $u_0(i Delta x)=c > 0$ for some $i$ and $u_0(j Delta x)=0$ $forall j != i$. Then:
      $
        u_i^1 = (1 - 2 lambda) u_i^0 + lambda (u_(i+1)^0 + u_(i-1)^0) = (1 - 2 lambda) c < 0
      $
      Now, the scheme is consistent and by the maximum principle for the continuous problem we have that $0 <= v(x,t)<=c$, where $v$ is the solution of the continuous problem. Then, by the Lax equivalence theorem, we cannot have convergence of the scheme, and therefore it is not stable. Thus, the stability condition for the explicit Euler scheme in the $L^oo$ norm is $lambda <= 1/2$.
    ]
]<ex:1>
#exercice[
  We now want to perform the same analysis of the $L^oo$ stability but now for the implicit Euler scheme. We want to show that
  $
    L_h (u) = (u_j^(n+1)-u_j^n)/(Delta t) - kappa (u_(j+1)^(n+1) - 2 u_j^(n+1) + u_(j-1)^(n+1))/(Delta x^2) = 0
  $ <eq:implicit>
  satisfies the discrete maximum principle, this time with no condition on $Delta x$ and $Delta t$. We consider Dirichlet boundary conditions, i.e. that formula @eq:implicit holds for $j in [2,N-1]$, and we impose $u_1^n = u_N^n = 0$ for all $n$. Show, by relating $u_j^(n+1)$ to $u_j^n$, that if $m <= 0 <= M$ are two constants such that
  $
    m <= u_j^0 <= M, quad forall j in [1,N]
  $
  then
  $
    m <= u_j^n <= M, quad forall j in [1,N] and n >= 0
  $
]
#solution[
  Doing a similar analysis as in @ex:1, we can write the implicit scheme as
  $
    bold(A)bold(u)^(n+1) = bold(u)^n
  $
  with 
  $
    bold(A) = mat(
      1 + 2 lambda, -lambda, 0,dots.c,0;
      -lambda, 1 + 2 lambda, -lambda, dots.down,dots.v;
      0, dots.down, dots.down, dots.down, 0;
      dots.v, dots.down, -lambda, 1 + 2 lambda, -lambda;
      0, dots.c, 0, -lambda, 1 + 2 lambda;
    ) = bold(I) + lambda bold(B)
  $
  where $bold(B)=mat(
    2, -1, 0, dots.c, 0;
    -1, 2, -1, dots.down, dots.v;
    0, dots.down, dots.down, dots.down, 0;
    dots.v, dots.down, -1,2,-1;
    0, dots.c, 0, -1, 2;
  )$. It suffices that the eigenvalues of $bold(A)^(-1)$ are positive and less than 1 to have the discrete maximum principle. It can be seen (check @tridial) that the eigenvalues of the tridiagonal matrix $bold(B)$ are 
  $
    lambda_k = 2 - 2 cos((k pi)/(N+1)), quad k=1,...,N
  $
  which implies that the eigenvalues of $bold(A)^(-1)$ are
  $
    mu_k = 1/(1 + lambda lambda_k) = 1/(1 + lambda (2 - 2 cos((k pi)/(N+1)))), quad k=1,...,N
  $
  which are positive and less or equal to 1. Thus, the implicit scheme satisfies the discrete maximum principle.
]
#exercice[
  We want to demosntrate Lax equivalence theorem for the heat equation using a linear finite difference scheme, with two levels. We introduce the solution $u(x,t)$ (supposed to be regular) of the heat equation @eq:heat, and $u_j^n$ the numerical approached solution using finite differences with $u_j^0=u_0(x_j)$. We want to show that if the scheme is consistent and stable for a given norm $||dot.c||$, then the scheme is convergent in the sense that the error vector $e_j^n = u_j^n - u(x_j,t_n)$ satisfies
  $
    lim_(Delta x, Delta t -> 0) lr((sup_(n Delta t <= T) ||e^n||)) = 0 quad forall T > 0
  $
  We will suppose the existence and uniqueness of the solution to the heat equation and use the definition of stability @eq:stability.
]
#solution[
  From the hypothesis we can write our problem as:
  $
    bold(A)_1 bold(u)^(n+1) = bold(A)_2 bold(u)^n
  $
  where we assume that the matrix $bold(A)_1$ is invertible (otherwise the problem is ill-posed or it has no uniqueness of solutions). Without loss of generality we may assume $norm(bold(A)_1^(-1))<= C Delta t$ (we can assume this either to $bold(A)_1^(-1)$ or to $bold(A)_2$). Let $v$ be the exact solution of the continuous problem and $bold(v)$ the vector composed of the values of $v$ at the grid points. We define the truncation error at time $n$ as $bold(T)^n = bold(A)_1 bold(v)^(n+1) - bold(A)_2 bold(v)^n$. We have that
  $
    bold(u)^(n+1) &= bold(A)_1^(-1) bold(A)_2 bold(u)^n\
    bold(v)^(n+1) &=  bold(A)_1^(-1) lr((bold(A)_2 bold(v)^n + bold(T)^n))
  $
  Let $bold(B) := bold(A)_1^(-1) bold(A)_2$. We then have:
  $
    bold(u)^(n+1) - bold(v)^(n+1) &= bold(B) (bold(u)^n - bold(v)^n) - bold(A)_1^(-1) bold(T)^n \
    &<= bold(B) (bold(u)^n - bold(v)^n) - bold(A)_1^(-1) bold(T)^n \
    &<= bold(B)^2 (bold(u)^(n-1) - bold(v)^(n-1)) - bold(B) bold(A)_1^(-1) bold(T)^(n-1) - bold(A)_1^(-1) bold(T)^n \
    &<= bold(B)^(n+1) (bold(u)^0 - bold(v)^0) - sum_(k=0)^n bold(B)^k bold(A)_1^(-1) bold(T)^(n-k) \
    &= - sum_(k=0)^n bold(B)^k bold(A)_1^(-1) bold(T)^(n-k)
  $<eq:convergence>
  If we assume that the scheme is consistent of order $p$ in space and $q$ in time we have that $||bold(T^n)|| = Order((Delta t)^q + (Delta x)^p)$. Now, it can be seen that stability condition is equivalent to having $norm(bold(B)^k)<=K$ $forall k in NN$ and $K$ given in @eq:stability. Indeed, if we have @eq:stability, then:
  $
    bold(u)^k = bold(B) bold(u)^(k -1) = ... = bold(B)^k bold(u)^0 ==> norm(bold(B)^k bold(u)^0) = norm(bold(u)^k) <= K norm(bold(u)^0)
  $
  which from the definition of matrix norm implies that $norm(bold(B)^k)<=C$, and this is valid $forall k$ such that $k Delta t<=T$. Now suppose that we have $norm(bold(B)^k)<=K$, then:
  $
    norm(bold(u)^k)=norm(bold(B)^k bold(u)^0)<=K norm(bold(u)^0)
  $
  Thus, the scheme is stable. With that in mind, continuing from @eq:convergence, we have that
  $
    ||bold(u)^(n+1) - bold(v)^(n+1)|| &<= sum_(k=0)^n norm(bold(B)^k) norm(bold(A)_1^(-1)) norm(bold(T)^(n-k))\
    & <= K C Delta t sum_(k=0)^n norm(bold(T)^(n-k)) \
    & =  K C (n + 1) Delta t Order((Delta t)^q + (Delta x)^p) \
    & <= K C (T + Delta t) Order((Delta t)^q + (Delta x)^p)
  $ 
  which goes to zero as $Delta t$ and $Delta x$ go to zero at the same order as the consistency of the scheme.
]
#exercice(numbering: none, title:"Problem")[
  We now consider the following one dimensional problem
  $
    pdv(u, x, 2) =-f, quad forall x in (0,1)
  $ <eq:poisson>
  where $f in cal(C)^1([0,1], RR)$, with homogeneous Dirichlet boundary conditions $u(0) = u(1) = 0$. We now introduce a Finite Volume discretization of the $[0,1]$ interval:
  $
    (Omega_j)_(j=1,...,N) quad "defined by" quad Omega_j = (x_(j-1/2), x_(j+1/2)), 
  $
  with 
  $
    x_(1/2) = 0 < x_(3/2) < ... < x_(j-1/2) < x_(j+1/2) < ... < x_(N-1/2) < x_(N+1/2) = 1
  $
  Finally, we consider the dual grid composed of $N$ points $(x_j)_(j=1,...,N)$ such that
  $
    x_(1/2) = 0< x_1 < x_(3/2) < ... < x_(j-1/2) < x_j < x_(j+1/2) < ... < x_(N) < x_(N+1/2) = 1
  $
  We introduce $h_(j+1/2) = x_(j+1) - x_j$ and $h_j = x_(j+1/2) - x_(j-1/2)$, where for simplicity we set $x_0 =0$ and $x_(N+1) = 1$. 
  1. By integrating @eq:poisson over $Omega_j$, write the discrete scheme in the finite volume sense for the unknows $u_j$. We will define $f_j$ for the discrete right-hand side of @eq:poisson and we will introduce $tilde(F)_(j+1/2)$ to identify the numerical approximate fluxes in $x_(j+1/2)$.
    #solution[
      Recall that we define $u_j$ as the average of $u$ over $Omega_j$, i.e. $u_j = 1/h_j integral_(Omega_j) u(x) dd(x)$. Integrating @eq:poisson over $Omega_j$ we get:
      $
        integral_(Omega_j) pdv(u, x, 2) dd(x) &= -integral_(Omega_j) f dd(x) \
        pdv(u, x)(x_(j+1/2)) - pdv(u, x)(x_(j-1/2)) &= -overline(f)_j h_j
      $ <eq:fvpoisson>
      where $overline(f)_j = 1/h_j integral_(Omega_j) f dd(x)$.
      Using centered finite differences for the derivative, we get:
      $
        (u_(j+1) - u_j)/h_(j+1/2) - (u_j - u_(j-1))/h_(j-1/2) = -overline(f)_j h_j
      $
    ]
  2. Under which hypothesis on the mesh is the scheme consistent? (Important: we will not make this hypothesis below)
    #solution[
      From one side, $overline(f)_j$ tends to $f(x_j)$ as $h_j$ tends to zero. Now, let $h_j^+:= x_(j+1/2) - x_j$, $h_j^-:= x_j - x_(j-1/2)$, $h_(j+1/2)^+:=x_(j+1)-x_(j+1/2)$ and $h_(j+1/2)^-:=x_(j+1/2)-x_(j)$. Then, $h_j = h_j^+ + h_j^-$ and $h_(j+1/2) = h_(j+1/2)^+ + h_(j+1/2)^-$. Taylor-expanding arround $pdv(u,x)(x_j)$ we have:
      $
        pdv(u,x)(x_(j+1/2)) &= pdv(u,x)(x_j) + pdv(u,x,2)(x_j)h_j^+ + pdv(u,x,3)(x_j)(h_j^+)^2/2 + Order((h_j^+)^2) \
        pdv(u,x)(x_(j-1/2)) &= pdv(u,x)(x_j) - pdv(u,x,2)(x_j)h_j^- + pdv(u,x,3)(x_j)(h_j^-)^2/2 + Order((h_j^-)^3)
      $ 
      This implies:
      $
        pdv(u,x,2)(x_j) &= (pdv(u, x)(x_(j+1/2)) - pdv(u, x)(x_(j-1/2)))/h_j + 1/h_j 1/2 pdv(u,x,3)(x_j) (h_j^2 - 2 h_j h_j^-) + Order((h_j^+)^2) + Order((h_j^-)^2)\
         &= (pdv(u, x)(x_(j+1/2)) - pdv(u, x)(x_(j-1/2)))/h_j + 1/2 pdv(u,x,3)(x_j) (h_j - 2 h_j^-) + Order((h_j^+)^2) + Order((h_j^-)^2)\
        &= 1/h_j lr([(u_(j+1) - u_j)/h_(j+1/2) -1/2 pdv(u,x,2)(x_(j+1/2))(h_(j+1/2)^+-h_(j+1/2)^-) + Order((h_(j+1/2))^2)- \
        & #h(0.9cm) - (u_j - u_(j-1))/h_(j-1/2)+ 1/2 pdv(u,x,2)(x_(j-1/2))(h_(j-1/2)^+-h_(j-1/2)^-)+ Order((h_(j-1/2))^2)]) + Order(h_j)
      $
      Thus, in order to have a consistent scheme we need to have: 
      $
      h_(j+1/2)^+=h_(j+1/2)^-, quad forall j=0,...,N
      $
    ]
  3. We now want to prove uniqueness of the solution for the discrete system obtained in question 1. We note that this is equivalent to showing that the only solution to this problem with no right-hand side is uniformly zero $u_j=0$ $forall j$. To this end, we suggest to multiply the expression obtained at question 1 by $u_j$ in order to introduce squared quantities.
    #solution[
      We have the scheme:
      $
        (u_(j+1) - u_j)/h_(j+1/2) - (u_j - u_(j-1))/h_(j-1/2) = 0, quad forall j=1,...,N
      $
      Multiplying by $u_j$ we get:
      $
        (u_j u_(j+1)) / h_(j+1/2) - u_j ^ 2 lr((1/h_(j+1/2)+1/h_(j-1/2))) + (u_(j-1)u_j) / h_(j-1/2) = 0, quad forall j=1,...,N
      $
      Summing over $j$, and using the fact that $u_0 = u_(N+1) = 0$, we get:
      $
        0&=sum_(j=1)^N lr([(u_j u_(j+1)) / h_(j+1/2) - u_j ^ 2 lr((1/h_(j+1/2)+1/h_(j-1/2))) + (u_(j-1)u_j) / h_(j-1/2)])\
        &= sum_(j=1)^(N-1) (u_j u_(j+1)) / h_(j+1/2) - u_N ^ 2/h_(N+1/2) - sum_(j=1)^(N-1) u_j ^ 2/h_(j+1/2)- u_1 ^ 2/h_(1/2)- sum_(j=2)^N u_j^2/h_(j-1/2) + sum_(j=2)^N (u_(j-1)u_j) / h_(j-1/2)\
        &= - u_N ^ 2/h_(N+1/2) - u_1 ^ 2/h_(1/2) + 1/h_(j+1/2) sum_(j=1)^(N-1) lr((2 u_j u_(j+1) - u_j ^ 2 - u_(j+1) ^ 2))\
        &= - u_N ^ 2/h_(N+1/2) - u_1 ^ 2/h_(1/2) - 1/h_(j+1/2) sum_(j=1)^(N-1) lr((u_(j+1) - u_j)^2)
      $
      which implies that all the terms are zero, and thus $u_j = 0$ $forall j$.
    ]
  4. Let $u in cal(C)^2([0,1], RR)$ be a solution of @eq:poisson. Introducing 
    $
      cal(F)_(j+1/2) = -pdv(u,x)(x_(j+1/2))
    $
    the exact flux in $x_(j+1/2)$ and $h=max_j(h_j)$ show the "consistancy of fluxes", in the sense that:
    $
    |cal(F)_(j+1/2) - tilde(F)_(j+1/2)| <= C_1 h quad "with" C_1 > 0
    $
    Under which hypothesis on the mesh is the approximation 2nd order? That is to say, satisfying
    $
      |cal(F)_(j+1/2) - tilde(F)_(j+1/2)| <= C_2 h^2 quad "with" C_2 > 0
    $
    #solution[
      Using @eq:fvpoisson the numerical fluxes are 
      $
        tilde(F)_(j+1/2) = -(u_(j+1) - u_j)/h_(j+1/2), quad forall j=1,...,N
      $
      and we define $F_(j+1/2) := -(u(x_(j+1)) - u(x_j))/h_(j+1/2)$. Recall the definitions of $h_(j+1/2)^+$ and $h_(j+1/2)^-$. Then:
      $
        u(x_(j+1)) &= u(x_(j+1/2)) + pdv(u,x)(x_(j+1/2))h_(j+1/2)^+ + 1/2 pdv(u,x,2)(xi_(j+1/2))( h_(j+1/2)^+)^2\
        u(x_j) &= u(x_(j+1/2)) -  pdv(u,x)(x_(j+1/2))h_(j+1/2)^- + 1/2 pdv(u,x,2)(eta_(j+1/2)) (h_(j+1/2)^-)^2
      $
      for some $xi_(j+1/2) in (x_(j+1/2), x_(j+1))$ and $eta_(j+1/2) in (x_j, x_(j+1/2))$. Substracting the equations and dividing by $h_(j+1/2)$ we get:
      $
        lr(|cal(F)_(j+1/2) - F_(j+1/2)|)&= lr(|pdv(u,x)(x_(j+1/2)) -(u(x_(j+1)) - u(x_j))/h_(j+1/2)|) \
        &= 1/h_(j+1/2)lr(|C_(11)(h_(j+1/2)^+)^2 - C_12 (h_(j+1/2)^-)^2|)\
        &<=1/h_(j+1/2)lr((tilde(C)_11 [(h_(j+1/2)^+)^2 - (h_(j+1/2)^-)^2]+ tilde(C)_12 (h_(j+1/2)^-)^2))\
        &<= tilde(C)_11 (h_(j+1/2)^+ - h_(j+1/2)^-) + tilde(C)_12 h_(j+1/2)^-\
        &<= C_1 h
      $
      for some $C_1,tilde(C)_11,tilde(C)_12 > 0$. If we want to prove 2nd order consistency, we need $u in cal(C)^3([0,1],RR)$. In that case we have:
      $
        u(x_(j+1)) &= u(x_(j+1/2))  +  pdv(u,x)(x_(j+1/2))h_(j+1/2)^+ + 1/2 pdv(u,x,2)(x_(j+1/2))( h_(j+1/2)^+)^2  + Order((h_(j+1/2)^+)^3)\
        u(x_j) &= u(x_(j+1/2))  -  pdv(u,x)(x_(j+1/2))h_(j+1/2)^- + 1/2 pdv(u,x,2)(x_(j+1/2)) (h_(j+1/2)^-)^2  + Order((h_(j+1/2)^-)^3)
      $
      Substracting the equations, we notice that to cancel out the second order term we need to have $h_(j+1/2)^+ = h_(j+1/2)^-$. The terms for the 3rd order coefficients follow as in the previous case.
    ]
  5. We now want to show convergence of the scheme (even when consistancy is not met). We suppose that the solution of @eq:poisson is in $cal(C)^2([0,1],RR)$ and introduce the error $e_j = u(x_j)-u_j$ for $j=1,...,N$, $e_0=e_(N+1)=0$. Using the consistancy of fluxes and Cauchy-Schwarz inequality, show that there exists $C_3>=0$ such that
    $
      sum_(j=0)^N ((e_(j+1)-e_j)^2)/h_(j+1/2) <= C_3^2 h^2
    $
    #solution[
      From the previous work we have the following equalities:
      $
        cal(F)_(j+1/2) - cal(F)_(j-1/2) &= overline(f)_j h_j\
        tilde(F)_(j+1/2) - tilde(F)_(j-1/2) &= overline(f)_j h_j\
        cal(F)_(j+1/2) - F_(j+1/2) &= T_(j+1/2) quad "with" |T_(j+1/2)| <= C_1 h
      $ <eq:Fs>
      Now we have that:
      $
        F_(j+1/2) - tilde(F)_(j+1/2) &= - (e_(j+1) - e_j)/h_(j+1/2) 
      $ <eq:errorF>
      Thus, putting together @eq:Fs and @eq:errorF we get:
      $
        (e_(j+1) - e_j)/h_(j+1/2) - (e_(j) - e_(j-1))/h_(j-1/2) = -T_(j+1/2)+T_(j-1/2), quad forall j=1,...,N
      $
      Multiplying this last expression by $e_j$ and summing over all $j$ we get:
      $
        sum_(j=1)^N ((e_(j+1)-e_j)e_j)/h_(j+1/2) - sum_(j=1)^N ((e_(j)-e_(j-1))e_j)/h_(j-1/2) &= -sum_(j=1)^N T_(j+1/2)e_j + sum_(j=1)^N T_(j-1/2)e_j\
        (e_(N+1) - e_N)e_N/h_(N+1/2) - sum_(j=1)^(N-1)((e_(j+1)-e_j)e_j)/h_(j+1/2) - &(e_1 - e_0)e_1/h_(1/2) - sum_(j=1)^(N-1) ((e_(j+1)-e_j)e_j)/h_(j+1/2) =\ 
        &= -T_(N+1/2)e_N -sum_(j=1)^N T_(j+1/2)e_j + T_(1/2) e_1+ sum_(j=1)^(N-1) T_(j+1/2)e_(j+1) \
        -e_N^2/(h_(N+1/2)) - e_1^2/(h_(1/2)) - sum_(j=1)^(N-1) ((e_(j+1)-e_j)^2)/h_(j+1/2) &= -T_(N+1/2)e_N + T_(1/2) e_1 + sum_(j=1)^(N-1) T_(j+1/2)(e_(j+1)-e_j)\
        sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2) &= - sum_(j=0)^(N) T_(j+1/2)(e_(j+1)-e_j)\
      $
      Finally taking absolutes values an using Cauchy-Schwarz inequality we get:
      $
        sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2) &<= C_1 h sum_(j=0)^(N) (|e_(j+1)-e_j|)/ sqrt(h_(j+1/2)) sqrt(h_(j+1/2))\
        &<= C_1 h sqrt(sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2)) sqrt(sum_(j=0)^(N) h_(j+1/2))\
        &= C_1 h sqrt(sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2))
      $
      where we have used that $sum_(j=0)^(N) h_(j+1/2) = 1$.
      Thus:
      $
        sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2) &<= C_1^2 h^2
      $
    ]
  6. Noting that 
    $
    |e_j| = lr(|sum_(k=1)^j e_k -e_(k-1)|)
    $
    conclude that
    $
      |e_j| <= C_3 h quad "for" j=1,...,N
    $
    which implies convergence of the finite volume scheme for @eq:poisson.
    #solution[
        Fix $i in {1,...,N}$. We have that:
        $
          |e_i| = lr(|sum_(j=1)^i e_j -e_(j-1)|) <= sum_(j=1)^i (|e_(j+1)-e_j|)/sqrt(h_(j+1/2)) sqrt(h_(j+1/2)) <= sqrt(sum_(j=1)^i ((e_(j+1)-e_j)^2)/h_(j+1/2)) sqrt(sum_(j=1)^i h_(j+1/2)) <= C_1 h
        $
        again by the Cauchy-Schwarz inequality, $sum_(j=1)^i h_(j+1/2) <= 1$ and the fact that:
        $
    sum_(j=1)^i ((e_(j+1)-e_j)^2)/h_(j+1/2)<= sum_(j=0)^(N) ((e_(j+1)-e_j)^2)/h_(j+1/2) <= C_1^2 h^2
        $
    ]
]
#bibliography("refs.bib",style: "ieee")
