%==============================================================================
% Daniel J. Greenhoe
% Experiments with odd polynomials R(1/2-y)
%
% References:
%    C. Sidney Burrus, Ramesh A. Gopinath, and Haitao Guo
%    Introduction to Wavelets and Wavelet Transforms; A Primer
%    Prentice Hall, 1998
%    ISBN 0-13-489600-9
%==============================================================================

%======================================
% Initialization
%======================================
clear;

%======================================
% Function Declarations
%======================================

%----------------------------------------------------------------------------
% \brief Downsample x(n) by L
% \params[in] x: data sequence
% \params[in] L: Downsample factor
% \returns Downsampled sequence y
%----------------------------------------------------------------------------
function y = downSample( x, L )
  n = length(x);
  m = ceil(n/L);
  y = 0*x(1:m);
  for i=0:(m-1)
    y(i+1) = x(i*L+1);
  endfor
endfunction

%----------------------------------------------------------------------------
% \brief Upsample x(n) by L
% \params[in] x: data sequence
% \params[in] L: Upsample factor
% \returns Upsampled sequence y
%----------------------------------------------------------------------------
function y = upSample( x, L )
  n = length(x);
  m = n * L;
  y = zeros(1,m);
  for i=0:(n -1)
    y(i*L+1) = x(i+1);
  endfor
endfunction

%----------------------------------------------------------------------------
% \brief Dilation Operation
% \details The dilation operator D on a function f(x) is defined as
%          D f(x) = sqrt(2) f(2x) .
%          The sqrt(2) factor makes D unitary.
% \params[in] x: data sequence
% \returns Dilated sequence y
%----------------------------------------------------------------------------
function y = Dilation(x)
  y = sqrt(2) * downSample(x, 2);
endfunction

%----------------------------------------------------------------------------
% \brief Make a sequence h(n) satisfy the Admissibility Condition
% \details A sequence h(n) satisfies the Admissibility Condition
%          if the sum of the elements of h(n) equals sqrt(2)
% \params[in] h: sequence
% \returns
%----------------------------------------------------------------------------
function hh = Admissibility(h)
  sumh = sum(h);
  hh = h;
  if( sumh==0 )
    printf("Admissibility(h) ERROR: sum(h)=0\n");
  else
    hh = (sqrt(2) / sumh ) * hh;
  end
endfunction

%----------------------------------------------------------------------------
% \brief Perform Dilation Equation operation
% \params[in]  x: data sequence
% \params[in]  h: scaling coefficients h(n)
% \returns Dilated sequence y
%----------------------------------------------------------------------------
function y = DilationEqu( x, h )
  span    =  length(h)-1;
  density = length(x) / span;          % number of x elements per h element except last
  N       = span * density;            % target length of y
  hu      = upSample(h,density);       % upsample h(n) to match phi(x) density
  xh      = conv(hu, x);               % phi_{k+1}(x) = SUM h(n) phi_k(2x-n)
  Dxh     = Dilation(xh);              % Dilation operation
  y       = Dxh(1:N);                  % truncate to length N
endfunction

%----------------------------------------------------------------------------
% \brief Estimate the integral of x
% \params[in]  x: data sequence
% \params[in]  density: number of elements per length 1 on x-axis
% \returns Estimated integral
%----------------------------------------------------------------------------
function INTx = Integrate( x, density )
  dx   = 1 / density;                  % width of Riemann integration bin
  INTx = sum( x ) * dx;                % estimated integral of phi(x)
endfunction

%----------------------------------------------------------------------------
% \brief Generate scaling function phi(x) from scaling coefficients h(n)
% \params[in]  h:          scaling coefficients h(n)
% \params[in]  iterations: number of iterations of Dilation Equation
% \params[in]  density:    number phi(x) samples per h(n) sample
% \returns     phi:        scaling function phi(x)
% \detailed   Cascade Algorithm: Strang and Nguyen (1996)
% reference: Burrus page 67
% reference: Rao
%----------------------------------------------------------------------------
function phi = gen_phi(h, iterations, density, verbose=1)
  h = Admissibility(h);                % Ensure h(n) satisfies Admiss. Cond.
  span =  length(h)-1                  % support of phi(x) = span of h(n)
  N    = span * density;               % target length of psi(x)
  phi = ones(1,density*span)/span;     % initial phi_0 is box with area=1
  for i = 0:(iterations-1)             % iterate Dilation Equation k times
    p0  = phi;                         % store phi
    phi = DilationEqu( phi, h );
    errorVect = phi - p0;
    cost = errorVect * errorVect';
   %if(verbose) printf("%f ", cost); end
  endfor
  if( verbose==1 )
    Iphi   = Integrate( phi, density );% estimated integral of psi(x)
    sumSq  = sum(h)^2;
    minPhi = min(phi);
    maxPhi = max(phi);
    printf("psi = gen_psi( phi, g, density=%d, verbose=1 )\n", density );
    printf("  * sum^2(h) = %12.8f", sumSq );   Test_equal( sumSq  , 2 );
    printf("  * min(phi) = %12.8f", minPhi );  Test_lt(    minPhi , 0 );
    printf("  * max(phi) = %12.8f", maxPhi );  Test_gt(    maxPhi , 0 );
    printf("  * INTphi   = %12.8f", Iphi);     Test_equal( Iphi   , 1 );
    printf("  * cost     = %12.8f", cost);     Test_equal( cost   , 0 );
%    printf("\nsum^2(h)=%f  min(phi)=%f  max(phi)=%f  INTphi=%f\n", sum(h)^2, min(phi), max(phi), Iphi);
  end
endfunction

%----------------------------------------------------------------------------
% \brief Generate wavelet function psi(x) from scaling function phi(x)
%        and wavelet filter coefficients g(n)
% \params[in] phi = scaling function phi(x)
% \params[in] g = wavelet filter coefficients g(n)
% \params[in] density = number phi(t) samples per g(n) sample
% \returns wavelet function psi(x)
% \cite Burrus page 15
%----------------------------------------------------------------------------
function psi = gen_psi(phi, g, density, verbose=1)
  span  = length(g)-1;                 % support of psi(x) = span of g(n)
  N     = span * density;              % target length of psi(x)
  psi   = DilationEqu( phi, g );       % psi(x) is a linear combination of phi(x)
  if( verbose )
    Ipsi   = Integrate( psi, density );% estimated integral of psi(x)
    sumSq  = sum(g)^2;
    minPsi = min(psi);
    maxPsi = max(psi);
    printf("psi = gen_psi( phi, g, density=%d, verbose=1 )\n", density );
    printf("  * sum^2(g) = %12.8f", sumSq );   Test_zero( sumSq );
    printf("  * min(psi) = %12.8f", minPsi );  Test_lt( minPsi, 0 );
    printf("  * max(psi) = %12.8f", maxPsi );  Test_gt( maxPsi, 0 );
    printf("  * INTpsi   = %12.8f", Ipsi);     Test_zero( Ipsi );
  end
endfunction

%----------------------------------------------------------------------------
% \brief Generate wavelet coefficients g(n) from scaling coefficients h(n)
% \params[in] h: scaling coefficients
% \returns Wavelet coefficients g
% \cite Mallat page 238
% \cite Burrus page  15
% \cite Burrus page  79
%----------------------------------------------------------------------------
function g = h2g_coefs(h)
  g = h;
  N = length(h);
  for n=0:(N-1)
     g(n +1) = (-1)**(n)*h(N-1-n +1);
  endfor
endfunction

%----------------------------------------------------------------------------
%function: Generate filter coefficients from wavelet/scalar coefficients
%   _
%   h(n) = h(-n) with shift by N where N = length of h(n)
%   _
%   g(n) = g(-n) with shift by N
%----------------------------------------------------------------------------
function hbar = h2hbar(h)
  hbar = zeros(size(h));               %allocate memory
  N    = length(h)     ;               %N = length of h(n)
  for n = 0:(N-1)                      %
     hbar(n +1) = h(N-1-n +1);         %hbar(n) = h(N-1-n)
  endfor
endfunction

%----------------------------------------------------------------------------
%                                               -z + 2 -1/z
% Convert polynomial in y to polynomial in y = -------------
%                                                    4
%                               /              \
%        |                     |  -z + 2 -1/z  |
%    P(y)|                  = P| ------------- |  = Q(z)Q(1/z)
%        |y=[(-z+2-1/z)/4]      \      4       /
%
% Input
% ------
%  P = [p_{n-1} ... p_2 p_1 p_0]
%    ==> P(y)=p_{n-1}y^{m-1}+...+p_2 y^2+p_1y+p_0
%
% Output
% ------
%  QQ = [q_{2n-2} ... q_2 q_1 q_0]
%
%                     q_{2n-2}z^{2n-2} +...+ q_2 z^2 + r_1 z + r_0
%     ==>  Q(z)Q(1/z)=---------------------------------------------
%                                       z^n
%
% Theory
% ------
% This conversion is useful when using polynomials in y = sin^2(w/2).
% Polynomials in sin^2(w/2) can be used to represent any even
% periodic funtion with period 2pi.
%   * sin^2(w/2) = (1/2)(1-cosw)
%
%                   1-cosw     1     e^{w}+e^{-w}     2 - z - 1/z   |
% y = sin^2(w/2) = -------- = --- - -------------- = -------------- |
%                      2       2          4                 4       |z=e^w
%
%----------------------------------------------------------------------------
function QQ = y2sin2(P)
  n=length(P);                         % number of terms in P(y). num 0s=n-1
  N=2*n-1;                             % number of terms in Q(z)Q(1/z)
  QQ= zeros(1,N);                      % init Q(z)Q(1/z) = 0+0z+...+0z^{N-1}
  q=1;                                 % init q(z)       = 1
  for k=0:n-1                          %
    QQ = [QQ,0](2:N+1);                % z[-z + 2 -1/z] = -z^2 + 2z - 1
    QQ = QQ + ...                      % Q(z)Q(1/z)
      P(n-k)*[zeros(1,N-2*k-1),q]/4^k; %   = SUM p1k_k*[(-z+2-z^-1)/4]^k
    q = conv(q,[-1 2 -1]);             % q(z) = q(z)[-z^2 + 2z - 1]
  endfor                               %
endfunction

%----------------------------------------------------------------------------
% Generate Daubechies class scaling coefficients
%  - Daubechies-p: minimum support (2p-1, R(y)=0) and minimum phase
%  - Symmlets-p:   minimum support (2p-1, R(y)=0) and quasi-linear phase
%  - Experimental: non-minimum support (>2p-1, R(y)!=0)
%
%  Input
%  ------
%  p = number of vanishing moments
%  R = [r_{m-1} ... r_2 r_1 r_0] ==> R(y)=r_{m-1}y^{m-1}+...+r_2 y^2+r_1y+r_0
%
%  Output
%  ------
%  n  = p + m
%
%  QQ = [q_{2n-1} ... q_1 q_0]
%
%                    q_{2n-1}z^{2n-1}+...+r_2 z^2+r_1z+r_0
%    ==>  Q(z)Q(1/z)=--------------------------------------
%                                  z^n
%
%  rQQ= [ r_1 r_2 ... r_{n-1}  r_n r_{n+1}...r_{n-2} ]
%        |<--roots inside -->|<--roots outside  --->|
%        |   unit circle     |   unit circle        |
%
%  A  = [a_p ... a_2 a_1 a_0]
%
%                       (z+1)^p
%    ==>  A(z) = sqrt(2)(---)   = a_pz^p + ... + a_2z^2 + a_1z + a_0
%                       ( 2 )
%
%  rA = [ r_1 r_2 ... r_p ] = roots of A(z)
%
%
%  Theory
%  ------
%  h(n) = scaling coefficients
%
%  H(z) = SUM h(n) z^{-n}
%          n
%       = A(z)Q(z)
%
%                (z+1)^p
%       = sqrt(2)(---)   Q(z)
%                ( 2 )
%
%  Finding Q(z) involves factoring a polynomial P(...):
%
%  P(y) = Pm(y) + y^p R(1/2-y)
%
%  In general, R is any polynomial that satisfies R(1/2-y)=-R(y-1/2)
%    (R is an odd polynomial about y=1/2)
%  For Daubechies-p wavelets and symmlets, R(y)=0
%    (giving Daubechies-p and symmlets minimum support)
%
%          p-1 ( p-1+k )
%  Pm(y) = SUM (       ) y^k
%          k=0 (   k   )
%
%  Q(z)Q(1/z) = P([2-z-1/z]/4)
%----------------------------------------------------------------------------
function [n,A,QQ,rA,rQQ] = gen_Dclass(p,R)
                                       % Compute A(z) = sqrt(2)[ (z+1)/2 ]^p
                                       % ------------------------------------
  A=1;                                 % A(z) = (z+1)^0 = 1
  for k=1:p                            %
    A=conv(A,[1 1]);                   % A(z) = (z+1)^k  k=1,2,3,...,p
  endfor                               %
  A = (sqrt(2)/2^p)*A;                 % A(z) = sqrt(2)[(z+1)/2]^p
  rA = sort(roots(A))';                % rA   = roots of A(z)
                                       
                                       % Pm(y)
                                       % ------------------------------------
  Pm = zeros(1,p);                     % init Pm(y) = 0y^{p-1}+...+0y^2+0y+0
  for k=0:p-1                          %
    Pm(p-k) = bincoeff(p-1+k,k);       % Pm(y) = SUM {p-1+k choose k} y^k
  endfor                               %
                                       
                                       % P( [2-z-1/z]/4 )
                                       % ------------------------------------
  m=length(R);                         %
  P=[zeros(1,m),Pm] + [R,zeros(1,p)];  % P(y) = Pm(y) + y^p R(y)
  while P(1)==0                        % remove leading zero coefficient terms
    n=length(P);                       %   n = number of terms in P(y)
    P=P(2:n);                          %   remove leading zero coefficient
  endwhile                             %
  n=length(P);                         % number of terms in P(y). num 0s=n-1
  QQ = y2sin2(P);                      % Q(z)Q(1/z) = P(y)|_y={(-z+2-1/z)/4}
                                       
                                       % Compute H(z)=A(z)Q(z)
                                       % ------------------------------------
  rQQ = sort(roots(QQ))';              % roots of Q(z)Q(1/z) (p-1 roots of Q(z), p-1 roots of Q(1/z))
  rQ  = rQQ(1:n-1);                    % roots of Q(z) (p-1 roots inside unit circle)
  Q   = poly(rQ);                      % convert roots into poly Q(z)
  Q   = sign(real(Q)).*abs(Q);         % eliminate any extraneous imag. components
  H   = conv(A,Q);                     % H(z) = A(z)Q(z)
  rH  = sort(roots(H))';               % roots of H(z)
  h   = Admissibility(H);              % Ensure h(n) satisfies Admiss. Cond.
endfunction

%----------------------------------------------------------------------------
%function: Generate Daubechies-p scaling coefficients
%
%  Input
%  ------
%  p:   number of vanishing moments
%
%  Output
%  ------
%  n  = p + m
%  h:   scaling coefficients
%  rQQ: roots of Q(z)Q(z^-1)
%                |   |_______ p-1 roots outside unit circle
%                |___________ p-1 roots inside  unit circle
%  rH:  roots of H(z)         p roots at z=-1 and p-1 roots of Q(z)
%
%
%  Theory
%  ------
%  h(n) = scaling coefficients
%
%  H(z) = SUM h(n) z^{-n}
%          n
%       = A(z)Q(z)
%
%                (z+1)^p
%       = sqrt(2)(---)   Q(z)
%                ( 2 )
%
%----------------------------------------------------------------------------
function [h,rQQ,rH] = gen_Dp(p)
  R = [0];                             % for Daubechies-p, R(y)=0
  [n,A,QQ,rA,rQQ] = gen_Dclass(p,R);   %
  rQ  = rQQ(1:p-1);                    % roots of Q(z) (p-1 roots inside unit circle)
  Q   = poly(rQ);                      % convert roots into poly Q(z)
  Q   = sign(real(Q)).*abs(Q);         % eliminate any extraneous imag. components
  H   = conv(A,Q);                     % H(z) = A(z)Q(z)
  rH  = sort(roots(H))';               % roots of H(z)
  h   = Admissibility(H);              % Ensure h(n) satisfies Admiss. Cond.
endfunction

%----------------------------------------------------------------------------
%function: Generate Symmlets-p scaling coefficients
%
%  Input
%  ------
%  p:   number of vanishing moments
%
%  Output
%  ------
%  h:   scaling coefficients
%  rQQ: roots of Q(z)Q(z^-1)
%                |   |_______ p-1 roots outside unit circle
%                |___________ p-1 roots inside  unit circle
%  rH:  roots of H(z)         p roots at z=-1 and p-1 roots of Q(z)
%
%----------------------------------------------------------------------------
function [h,rQQ,rH] = gen_Sp(p)
                                       % Initialization
                                       % ---------------------
  N = 1024;                            % number of data points for linear measure
  R = [0];                             % for Symmlets-p, R(y)=0
  [n,A,QQ,rA,rQQ]=gen_Dclass(p,R);     %
  Nquads    = floor(p/2);              % number of zero quads
  Nsol      = 2^Nquads;                % number of possible solutions
  rQQ_mat    = zeros(p-1,Nsol);        % allocate matrix of roots for all solutions
  mse       = zeros(1,  Nsol);         % allocate mean square error of each solution
  sqMag     = zeros(1,  Nsol);         % allocate square magnitude of each solution

                                       % Main loop: generate all root selections
                                       % ---------------------------------------
  for n=0:Nsol-1                       % loop once for each possible solution
    t = conj(rQQ(1:2:p-1));            % extract one representative from each quad
    for m=0:Nquads-1                   % take conjugate recipricols of roots
      if( mod(floor(n/2^m),2) )        % if bit k of n is 1,
        t(m+1) = conj(1/t(m+1));       %   then take conjugate recipricol of root
      endif                            %
    endfor                             %

    rNew          = zeros(p-1,1);      % initialize new root column vector
    rNew(1:2:p-1) = t;                 % rNew=[r1 (real), r2, r2*, r3, r3*,...]'
    rNew(2:2:p-1) = conj(t(mod(p+1,2)+1:Nquads));
    rMat(:,n+1)   = rNew;              % store new root vector in root matrix

                                       % Linear phase measurements
                                       % ------------------------------
    h = real(poly(rNew));              % convert new roots into a poly in z
    [H,w] = freqz(h,[1],N);            % compute the frequency response
    phase = unwrap(angle(H));          % unwrapped phase of h
    [pCoefs,pVals]=polyfit(w,phase,1); % find the 1st order poly p(w) that best fits phase(w)
    mse(n+1) = (pVals.yf-phase)'*(pVals.yf-phase)/N;  % measure the error
    sqMag(n+1) = (rNew'*rNew);         % measure the magnitude and store
  endfor
                                       % Find best set of roots
                                       % --------------------------
  [err1,i1] = min(mse);                % find the 1st location of minimum error
  t=mse;                               % copy mse to temporary vector t
  t(i1) = max(t);                      % get 1st minimum out of the way to find 2nd
  [err2,i2] = min(t);                  % find the 2nd location of minimum error
  t= max(mse)*ones(1,Nsol);            % set t to all large values
  if( sqMag(i1)<sqMag(i2) )            % decide which min to use based on smaller magnitude
    iMinErr = i1;                      %   if the magnitude of first is smaller
  else                                 %
    iMinErr = i2;                      %   if the magnitude of the second is smaller
  endif                                %
  mmse = mse(iMinErr);                 % minimum mean square error
  rNew = rMat(:,iMinErr);              % extract min phase min magnitude vector

                                       % Compute H(z)=A(z)Q(z)
                                       % ------------------------------------
  rQQ = [rNew; conj(1./rNew)]';        % vector w/ both retained and discarded roots
  rQ  = rQQ(1:p-1);                    % roots of Q(z) (p-1 roots inside unit circle)
  Q   = poly(rQ);                      % convert roots into poly Q(z)
  Q   = sign(real(Q)).*abs(Q);         % eliminate any extraneous imag. components
  H   = conv(A,Q);                     % H(z) = A(z)Q(z)
  rH  = sort(roots(H))';               % roots of H(z)
  h   = Admissibility(H);              % Ensure h(n) satisfies Admiss. Cond.
endfunction

%----------------------------------------------------------------------------
% Generate experimental coefficients with R(y)!=0
%
%  Input
%  ------
%  p:   number of vanishing moments
%  R = [r_{m-1} ... r_2 r_1 r_0] ==> R(y)=r_{m-1}y^{m-1}+...+r_2 y^2+r_1y+r_0
%
%  Output
%  ------
%  h:   scaling coefficients
%  rQQ: roots of Q(z)Q(z^-1)
%                |   |_______ p-1 roots outside unit circle
%                |___________ p-1 roots inside  unit circle
%  rH:  roots of H(z)         p roots at z=-1 and p-1 roots of Q(z)
%
%
%  Theory
%  ------
%  h(n) = scaling coefficients
%
%  H(z) = SUM h(n) z^{-n}
%          n
%       = A(z)Q(z)
%
%                (z+1)^p
%       = sqrt(2)(---)   Q(z)
%                ( 2 )
%
%----------------------------------------------------------------------------
function [n,h,rQQ,rH] = gen_Rp(p,R)
  [n,A,QQ,rA,rQQ] = gen_Dclass(p,R);   %
  rQ  = rQQ(1:n-1);                    % roots of Q(z) (n-1 roots inside unit circle)
  Q   = poly(rQ);                      % convert roots into poly Q(z)
  Q   = sign(real(Q)).*abs(Q);         % eliminate any extraneous imag. components
  H   = conv(A,Q);                     % H(z) = A(z)Q(z)
  rH  = sort(roots(H))';               % roots of H(z)
  h   = Admissibility(H);              % Ensure h(n) satisfies Admiss. Cond.
endfunction

#----------------------------------------------------------------------------
#function: Generate scaling coefficients h(n) from length-4 parameter alpha
#          Pollen wavelets
#  references:
#    http://www2.isye.gatech.edu/~brani/datasoft/DL.pdf  page 2
#    Burrus page  66
#----------------------------------------------------------------------------
function h = gen_pollen4(alpha)
  h =  (sqrt(2)/4) * [...
       1 + cos(alpha) - sin(alpha) ;
       1 + cos(alpha) + sin(alpha) ;
       1 - cos(alpha) + sin(alpha) ;
       1 - cos(alpha) - sin(alpha) ;
      ];
%  h =  (sqrt(2)/4) * [\  Burrus ?
%       1 - cos(alpha) + sin(alpha) ;
%       1 + cos(alpha) + sin(alpha) ;
%       1 + cos(alpha) - sin(alpha) ;
%       1 - cos(alpha) - sin(alpha) ;
%      ];
endfunction

%-----------------------------------
% Write data to file
%-----------------------------------
function data2file(h,g,rR,rH,filename,comment)
  data  = fopen(filename,"w");
  fprintf(data,'%%========================================================================\n');
  fprintf(data,'%% Daniel J. Greenhoe \n');
  fprintf(data,'%% file: %s \n',filename);
  fprintf(data,'%% %s \n',comment);
  fprintf(data,'%% %s \n',strftime('%Y %B %d %A %r (%Z)',localtime(time)));
  fprintf(data,'%%========================================================================\n\n');

  p=length(h)/2;
  fprintf(data ,'scaling coefficients\n' );
  fprintf(data ,'-------------------------\n' );
  for n=1:length(h)
    fprintf(data ,'h_{%02d} = %13.10f  \n',n-1,h(n) );
  endfor

  fprintf(data ,'\n' );
  fprintf(data ,'wavelet coefficients\n');
  fprintf(data ,'-------------------------\n' );
  for n=1:length(g)
    fprintf(data ,'g_{%02d} = %13.10f  \n',n-1,g(n) );
  endfor

  fprintf(data ,'\n' );
  fprintf(data ,'roots of Q(z)\n' );
  fprintf(data ,'-------------------------\n' );
  for n=1:p-1
    fprintf(data ,'r_{%02d} = %13.10f  %13.10fi \n',n-1, real(rR(n)), imag(rR(n)) );
  endfor

  fprintf(data ,'\n' );
  fprintf(data ,'roots of Q(z^-1)\n' );
  fprintf(data ,'-------------------------\n');
  for n=p:length(rR)
    fprintf(data ,'r_{%02d} = %13.10f  %13.10fi \n',n-1, real(rR(n)), imag(rR(n)) );
  endfor

  fprintf(data ,'\n');
  fprintf(data ,'roots of H(z)\n');
  fprintf(data ,'-------------------------\n');
  for n=1:length(rH)
    fprintf(data ,'r_{%02d} = %13.10f  %13.10fi \n',n-1, real(rH(n)), imag(rH(n)) );
    if(n==p-1) fprintf(data,'\n'); endif
  endfor

  fprintf(data ,'\n' );
  fprintf(data ,'LaTeX drawing commands   \n');
  fprintf(data ,'-------------------------\n');
  for n=1:p-1
    fprintf(data ,'  \\put( %15.10f, %15.10f ) {\\circle {15}} \n',real(rR(n))*100, imag(rR(n))*100 );
  endfor
  fprintf(data,"\n");
  for n=p:length(rR)
    fprintf(data ,'  \\put( %15.10f, %15.10f ) {\\circle*{15}} \n',real(rR(n))*100, imag(rR(n))*100 );
  endfor

  fflush(data);
  fclose(data);
  printf('output file name = \"%s\" \n',filename);
endfunction

%----------------------------------------------------------------------------
% Write data to plot file
% <filename> can be used by the command '\fileplot' in a TeX environment
% to generate a plot of (x,y).
% The command '\fileplot' is available in the PSTricks package of TeX.
% Reference: http://www.ctan.org/pkg/pstricks
%----------------------------------------------------------------------------
function data2plotfile(x,y,filename,comment)
  data  = fopen(filename,'w');
  fprintf(data,'%%========================================================================\n');
  fprintf(data,'%% Daniel J. Greenhoe \n');
  fprintf(data,'%% file: %s \n',filename);
  fprintf(data,'%% number of data points = %d\n',length(x));
  fprintf(data,'%% %s \n',comment);
  fprintf(data,'%% %s \n',strftime('%Y %B %d %A %r (%Z)',localtime(time)));
  fprintf(data,'%% The command \\fileplot{%s} may be used in a LaTeX environment for plotting data.\n',filename);
  fprintf(data,'%% \\fileplot is available in the LaTeX PSTricks package.\n');
  fprintf(data,'%% Reference: http://www.ctan.org/pkg/pstricks\n');
  fprintf(data,'%%========================================================================\n');

  fprintf(data ,'[\n');
  for n=1:length(x)
    fprintf(data ,'(%13.10f,%13.10f)\n',x(n),y(n));
  endfor
  fprintf(data,']\n');
  fprintf(data,'%%---end of file---' );
  fflush(data);
  fclose(data);
  printf("data written to \'%s\' \n",filename);
endfunction

%----------------------------------------------------------------------------
% Demonstration of Daubechies-p wavelets
%  iterations: number of iterations to use to generate phi(t) (e.g. 16)
%  N:          data size (e.g. 1024)
%----------------------------------------------------------------------------
function demo_Dp(p,N,iterations)
                                       % calculate coefficient sequences
                                       % ------------------------
  [h,rQQ,rH] = gen_Dp(p);              % Daubechies-p
  g   = h2g_coefs(h);                  % generate wavelet coefficients g(n)
  d   = round(N/length(h));            % density=round(samples/unit(1))
  phi = gen_phi(h,  iterations,d);     % generate phi(x) from h(n)
  psi = gen_psi(phi,g,d);              % generate psi(x) from g(n)
  M   = length(phi);

                                       % Generate output data files
                                       % ------------------------
  x = [0:length(phi)-1]*(length(h)-1)/(length(phi)-1);
  data2file(h,g,rQQ,rH,sprintf('d%d.dat',p),    sprintf('data for Daubechies-%d wavelets',p) );
  data2plotfile(x,phi,sprintf('d%d_phi.dat',p),sprintf('plot file for Daubechies-%d scaling function',p));
  data2plotfile(x,psi,sprintf('d%d_psi.dat',p),sprintf('plot file for Daubechies-%d wavelet function',p));
  plot(x,phi,x,psi);

endfunction

%----------------------------------------------------------------------------
% Demonstration of Symmlet-p wavelets
%  iterations: number of iterations to use to generate phi(t) (e.g. 16)
%  N:          data size (e.g. 1024)
%----------------------------------------------------------------------------
function demo_Symmlet_p( p, N, iterations, dataDump=0 )
  [h,rQQ,rH] = gen_Sp(p);              % Symmlet-p
  g       = h2g_coefs(h);              % generate wavelet coefficients g(n)
  span    = length(h) - 1;             % span of phi(x) and psi(x) = span(h)
  density = round( N / span );         % density=round(samples/unit(1))
  phi     = gen_phi(h, iterations, density); % generate phi(x) from h(n)
  psi     = gen_psi(phi,g,density);    % generate psi(x) from g(n)
  M       = length(phi);
  xAxis   = [0:M-1] * span/(M-1);      % x-axis vector from 0 to span
  if( dataDump )
    data2file(h,g,rQQ,rH, sprintf("s%d.dat",p),    sprintf("Symmlet-%d data file",p));
    data2plotfile(xAxis, phi,  sprintf('s%d_phi.dat',p),sprintf('plot file for Symmlet-%d scaling function',p));
    data2plotfile(xAxis, psi,  sprintf('s%d_psi.dat',p),sprintf('plot file for Symmlet-%d wavelet function',p));
  end
  plot(xAxis, phi, xAxis, psi );
endfunction

%----------------------------------------------------------------------------
% Demonstration of Daubechies-p class wavelet with R(y)!=0
%  iterations: number of iterations to use to generate phi(t) (e.g. 16)
%  N:          data size (e.g. 1024)
%----------------------------------------------------------------------------
function demo_Ry_p(p,N,iterations)
                                       % R(y)
                                       % ------------------------
  R=[  3 0 5 0 7 0 3 0 1 0];                        % R(y)
  ptitle=sprintf("p=%d      R(y)=3y^9 + 5y^7 + 7y^5 + 3y^3 + y      Daniel J. Greenhoe",p);% title

                                       % Generate scaling coefficients {h_n}
                                       % ------------------------
  [n,h,rQQ,rH]= gen_Rp(p,R);           % Ry(p)
  g   = h2g_coefs(h);                  % generate wavelet coefficients g(n)
  d   = round(N/length(h));            % density=round(samples/unit(1))
  phi = gen_phi(h,  iterations,d);     % generate phi(x) from h(n)
  psi = gen_psi(phi,g,d);              % generate psi(x) from g(n)

                                       % Generate output
                                       % -------------------------
  data2file(h,g,rQQ,rH, sprintf("R%d.dat",p), sprintf("R(y)-%d data file",p) );   pause(delay_sec);
endfunction

%----------------------------------------------------------------------------
% Demonstration of Pollen length 4 wavelets
%  iterations: number of iterations to use to generate phi(t) (e.g. 16)
%  N:          data size (e.g. 1024)
%  alpha:      parameter alpha
%  p4_phi.dat and p4_psi.dat can be used with LaTeX's \plotfile command.
%  \pltfile is available in the pstricks package
%   Reference: http://www.ctan.org/pkg/pstricks
%----------------------------------------------------------------------------
function demo_pollen4a(alpha,N,iterations)
                                       % calculate coefficient sequences
                                       % ------------------------
  h   = gen_pollen4(alpha);            % Pollen length-4
  g   = h2g_coefs(h);                  % generate wavelet coefficients g(n)
  d   = round(N/length(h));            % density=round(samples/unit(1))
  phi = gen_phi(h,  iterations,d);     % generate phi(x) from h(n)
  psi = gen_psi(phi,g,d);              % generate psi(x) from g(n)
  M   = length(phi);

                                       % Generate output data files
                                       % ------------------------
  x = [0:length(phi)-1]*(length(h)-1)/(length(phi)-1);
  %data2file(h,g,rQQ,rH,'p4.dat',   sprintf('data for Pollen-4 alpha=%13.10f wavelets',alpha));
  data2plotfile(x,phi,'p4_phi.dat',sprintf('plot file for Pollen-4 alpha=%13.10f scaling function',alpha));
  data2plotfile(x,psi,'p4_psi.dat',sprintf('plot file for Pollen-4 alpha=%13.10f wavelet function',alpha));
  plot(x,phi,x,psi);

endfunction

%----------------------------------------------------------------------------
% Demonstration of Pollen length 4 wavelets
%  iterations: number of iterations to use to generate phi(t) (e.g. 16)
%  N:          data size (e.g. 1024)
%  nalpha:     number of alpha values
%  a:          starting alpha value
%  b:          ending alpha value
%  p4_phi.dat and p4_psi.dat can be used with GNU Plot's 'splot' function
%  Reference: http://www.gnuplot.info
%----------------------------------------------------------------------------
function demo_pollen4(a,b,nalpha,N,iterations)
  filename1='p4_psi.dat';
  filename2='p4_phi.dat';
  data1  = fopen(filename1,'w');
  data2  = fopen(filename2,'w');
  fprintf(data1,'#========================================================================\n');
  fprintf(data1,'# Daniel J. Greenhoe \n');
  fprintf(data1,'# file: %s \n',filename1);
  fprintf(data1,'# %s \n',strftime('%Y %B %d %A %r (%Z)',localtime(time)));
  fprintf(data1,'# The command splot{%s} may be used in a GNU Plot environment for plotting data.\n',filename1);
  fprintf(data1,'# Reference: http://www.gnuplot.info\n');
  fprintf(data1,'#========================================================================\n');

  fprintf(data2,'#========================================================================\n');
  fprintf(data2,'# Daniel J. Greenhoe \n');
  fprintf(data2,'# file: %s \n',filename2);
  fprintf(data2,'# %s \n',strftime('%Y %B %d %A %r (%Z)',localtime(time)));
  fprintf(data2,'# The command splot{%s} may be used in a GNU Plot environment for plotting data.\n',filename2);
  fprintf(data2,'# Reference: http://www.gnuplot.info\n');
  fprintf(data2,'#========================================================================\n');

  for i = 0:nalpha-1
                                       % calculate coefficient sequences
                                       % ------------------------
    alpha = a + i/(nalpha-1)*(b-a)     % calculate alpha
    h   = gen_pollen4(alpha);          % Pollen length-4
    g   = h2g_coefs(h);                % generate wavelet coefficients g(n)
    d   = round(N/length(h));          % density=round(samples/unit(1))
    phi = gen_phi(h,  iterations,d);   % generate phi(x) from h(n)
    psi = gen_psi(phi,g,d);            % generate psi(x) from g(n)
    M   = length(phi);

                                       % Generate output data files
                                       % ------------------------
    x = [0:length(phi)-1]*(length(h)-1)/(length(phi)-1);
    y1 = psi;
    y2 = phi;

    for n=1:length(x)
      fprintf(data1 ,'%13.10f  %13.10f  %13.10f\n',alpha,x(n),y1(n));
    endfor
    for n=1:length(x)
      fprintf(data2 ,'%13.10f  %13.10f  %13.10f\n',alpha,x(n),y2(n));
    endfor
    fprintf(data1 ,'\n');
    fprintf(data2 ,'\n');

  endfor
  fprintf(data1,'##---end of file---' );
  fprintf(data2,'##---end of file---' );
  fflush(data1);
  fflush(data2);
  fclose(data1);
  fclose(data2);
  printf("data written to \'%s\' \n",filename1);
  printf("data written to \'%s\' \n",filename2);

endfunction

%======================================
% Test Utilities
%======================================
%----------------------------------------------------------------------------
% \brief Test if x less than y
% \returns 1 if TRUE, 0 otherwise
%----------------------------------------------------------------------------
function result = Test_lt( x, y, verbose=1 )
  if( x < y ) result = 1;
  else        result = 0;
  end
  if( verbose==1 )
    if( result == 1 ) printf("  ok");
    else              printf("  FAIL!!!");
    end
    printf("\n");
  end
endfunction

%----------------------------------------------------------------------------
% \brief Test if x greater than y
% \returns 1 if TRUE, 0 otherwise
%----------------------------------------------------------------------------
function result = Test_gt( x, y, verbose=1 )
  result = Test_lt( y, x, verbose );
endfunction

%----------------------------------------------------------------------------
% \brief Test if x almost equals y
% \returns 1 if TRUE, 0 otherwise
%----------------------------------------------------------------------------
function result = Test_equal( x, y, verbose=1 )
  result = Test_lt( abs(x-y), 1e-9, verbose );
endfunction

%----------------------------------------------------------------------------
% \brief Test if x is close to 0
% \returns 1 if TRUE, 0 otherwise
%----------------------------------------------------------------------------
function result = Test_zero( x, verbose=1 )
  result = Test_lt( abs(x), 1e-9, verbose );
endfunction

%======================================
% Main
%======================================
                                       % parameters
                                       %-------------------------------------
N = 1024;                              % number of data points
iterations = 100;                       % number of iterations

                                       % demos
                                       %-------------------------------------
%for p=1:12
%  demo_Dp(p,N,iterations);
%endfor
%demo_Dp(16,N,iterations);

%demo_Symmlet_p( 4,N,iterations);
%demo_Symmlet_p( 8,N,iterations);
%demo_Symmlet_p(12,N,iterations);
%demo_Symmlet_p(16,N,iterations);
%demo_Ry_p(3);
%demo_pollen4a(pi/6,16,iterations)%
%demo_pollen4a(pi/4,16,iterations)%
%demo_pollen4a(pi/2,16,iterations)%
%demo_pollen4a(pi,16,iterations)%
%demo_pollen4(0,pi,256,256,iterations)%
%demo_pollen4(0,pi,32,32,iterations)%
%function demo_pollen4(a,b,nalpha,N,iterations)

%h = [1 3 3 1]%
%phi = gen_phi(h,3,4)
%plot(phi)
%======================================
% End Processing
%======================================
p=4
  [h,rQQ,rH] = gen_Dp(p);              % Daubechies-p
  g   = h2g_coefs(h);                  % generate wavelet coefficients g(n)
  d   = round(N/length(h));            % density=round(samples/unit(1))
  phi = gen_phi(h,  iterations,d);     % generate phi(x) from h(n)
  psi = gen_psi(phi,g,d);              % generate psi(x) from g(n)
  M   = length(phi);
