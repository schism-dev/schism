!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Spherepack                            *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! PURPOSE
! -------
!     This package consists of programs which perform fast fourier
!     transforms for both complex and real periodic sequences and
!     certain other symmetric sequences that are listed below.
!
! USAGE
! -----
!     1.   rffti     initialize  rfftf and rfftb
!     2.   rfftf     forward transform of a real periodic sequence
!     3.   rfftb     backward transform of a real coefficient array
!
!     4.   ezffti    initialize ezfftf and ezfftb
!     5.   ezfftf    a simplified real periodic forward transform
!     6.   ezfftb    a simplified real periodic backward transform
!
!     7.   sinti     initialize sint
!     8.   sint      sine transform of a real odd sequence
!
!     9.   costi     initialize cost
!     10.  cost      cosine transform of a real even sequence
!
!     11.  sinqi     initialize sinqf and sinqb
!     12.  sinqf     forward sine transform with odd wave numbers
!     13.  sinqb     unnormalized inverse of sinqf
!
!     14.  cosqi     initialize cosqf and cosqb
!     15.  cosqf     forward cosine transform with odd wave numbers
!     16.  cosqb     unnormalized inverse of cosqf
!
!     17.  cffti     initialize cfftf and cfftb
!     18.  cfftf     forward transform of a complex periodic sequence
!     19.  cfftb     unnormalized inverse of cfftf
!
! SPECIAL CONDITIONS
! -------
!     before calling routines ezfftb and ezfftf for the first time, 
!     or before calling ezfftb and ezfftf with a different length, 
!     users must initialize by calling routine ezffti.
!
! HISTORY
! -------
!     * Developed at NCAR in boulder, colorado by paul n. swarztrauber
!       of the scientific computing division.  Released on NCAR's public
!       software libraries in January 1980.
!     * Modified may 29, 1985 to increase efficiency.
!     * Updated by Jon Lo Kim Lin in 2016 to incorporate features
!       of Fortran 2008
!
! **********************************************************************
!
!     subroutine rffti(n, wsave)
!
!     subroutine rffti initializes the array wsave which is used in
!     both rfftf and rfftb. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 2*n+15.
!             the same work array can be used for both rfftf and rfftb
!             as long as n remains unchanged. different wsave arrays
!             are required for different values of n. the contents of
!             wsave must not be changed between calls of rfftf or rfftb.
!
! **********************************************************************
!
!     subroutine rfftf(n, r, wsave)
!
!     subroutine rfftf computes the fourier coefficients of a real
!     perodic sequence (fourier analysis). the transform is defined
!     below at output parameter r.
!
!     input parameters
!
!     n       the length of the array r to be transformed.  the method
!             is most efficient when n is a product of small primes.
!             n may change so long as different work arrays are provided
!
!     r       a real array of length n which contains the sequence
!             to be transformed
!
!     wsave   a work array which must be dimensioned at least 2*n+15.
!             in the program that calls rfftf. the wsave array must be
!             initialized by calling subroutine rffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by rfftf and rfftb.
!
!
!     output parameters
!
!     r       r(1) = the sum from i=1 to i=n of r(i)
!
!             if n is even set l =n/2   , if n is odd set l =(n+1)/2
!
!               then for k = 2, ..., l
!
!                  r(2*k-2) = the sum from i = 1 to i = n of
!
!                       r(i)*cos((k-1)*(i-1)*2*pi/n)
!
!                  r(2*k-1) = the sum from i = 1 to i = n of
!
!                      -r(i)*sin((k-1)*(i-1)*2*pi/n)
!
!             if n is even
!
!                  r(n) = the sum from i = 1 to i = n of
!
!                       (-1)**(i-1)*r(i)
!
!      *****  note
!                  this transform is unnormalized since a call of rfftf
!                  followed by a call of rfftb will multiply the input
!                  sequence by n.
!
!     wsave   contains results which must not be destroyed between
!             calls of rfftf or rfftb.
!
!
! **********************************************************************
!
!     subroutine rfftb(n, r, wsave)
!
!     subroutine rfftb computes the real perodic sequence from its
!     fourier coefficients (fourier synthesis). the transform is defined
!     below at output parameter r.
!
!     input parameters
!
!     n       the length of the array r to be transformed.  the method
!             is most efficient when n is a product of small primes.
!             n may change so long as different work arrays are provided
!
!     r       a real array of length n which contains the sequence
!             to be transformed
!
!     wsave   a work array which must be dimensioned at least 2*n+15.
!             in the program that calls rfftb. the wsave array must be
!             initialized by calling subroutine rffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by rfftf and rfftb.
!
!
!     output parameters
!
!     r       for n even and for i = 1, ..., n
!
!                  r(i) = r(1)+(-1)**(i-1)*r(n)
!
!                       plus the sum from k=2 to k=n/2 of
!
!                        TWO * r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                       -TWO * r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!             for n odd and for i = 1, ..., n
!
!                  r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
!
!                       TWO * r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                      -TWO * r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!      *****  note
!                  this transform is unnormalized since a call of rfftf
!                  followed by a call of rfftb will multiply the input
!                  sequence by n.
!
!     wsave   contains results which must not be destroyed between
!             calls of rfftb or rfftf.
!
!
! **********************************************************************
!
!     subroutine ezffti(n, wsave)
!
!     subroutine ezffti initializes the array wsave which is used in
!     both ezfftf and ezfftb. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             the same work array can be used for both ezfftf and ezfftb
!             as long as n remains unchanged. different wsave arrays
!             are required for different values of n.
!
!
! **********************************************************************
!
!     subroutine ezfftf(n, r, azero, a, b, wsave)
!
!     subroutine ezfftf computes the fourier coefficients of a real
!     perodic sequence (fourier analysis). the transform is defined
!     below at output parameters azero, a and b. ezfftf is a simplified
!     but slower Version of rfftf.
!
!     input parameters
!
!     n       the length of the array r to be transformed.  the method
!             is must efficient when n is the product of small primes.
!
!     r       a real array of length n which contains the sequence
!             to be transformed. r is not destroyed.
!
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             in the program that calls ezfftf. the wsave array must be
!             initialized by calling subroutine ezffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by ezfftf and ezfftb.
!
!     output parameters
!
!     azero   the sum from i=1 to i=n of r(i)/n
!
!     a, b     for n even b(n/2)=0. and a(n/2) is the sum from i=1 to
!             i=n of (-1)**(i-1)*r(i)/n
!
!             for n even define kmax=n/2-1
!             for n odd  define kmax=(n-1)/2
!
!             then for  k=1, ..., kmax
!
!                  a(k) equals the sum from i=1 to i=n of
!
!                       2./n*r(i)*cos(k*(i-1)*2*pi/n)
!
!                  b(k) equals the sum from i=1 to i=n of
!
!                       2./n*r(i)*sin(k*(i-1)*2*pi/n)
!
!
! **********************************************************************
!
!     subroutine ezfftb(n, r, azero, a, b, wsave)
!
!     subroutine ezfftb computes a real perodic sequence from its
!     fourier coefficients (fourier synthesis). the transform is
!     defined below at output parameter r. ezfftb is a simplified
!     but slower Version of rfftb.
!
!     input parameters
!
!     n       the length of the output array r.  the method is most
!             efficient when n is the product of small primes.
!
!     azero   the constant fourier coefficient
!
!     a, b     arrays which contain the remaining fourier coefficients
!             these arrays are not destroyed.
!
!             the length of these arrays depends on whether n is even or
!             odd.
!
!             if n is even n/2    locations are required
!             if n is odd(n-1)/2 locations are required
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             in the program that calls ezfftb. the wsave array must be
!             initialized by calling subroutine ezffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by ezfftf and ezfftb.
!
!
!     output parameters
!
!     r       if n is even define kmax=n/2
!             if n is odd  define kmax=(n-1)/2
!
!             then for i=1, ..., n
!
!                  r(i)=azero plus the sum from k=1 to k=kmax of
!
!                  a(k)*cos(k*(i-1)*2*pi/n)+b(k)*sin(k*(i-1)*2*pi/n)
!
!     ********************* complex notation **************************
!
!             for j=1, ..., n
!
!             r(j) equals the sum from k=-kmax to k=kmax of
!
!                  c(k)*exp(i*k*(j-1)*2*pi/n)
!
!             where
!
!                  c(k) = .5*cmplx(a(k), -b(k))   for k=1, ..., kmax
!
!                  c(-k) = conjg(c(k))
!
!                  c(0) = azero
!
!                       and i=sqrt(-1)
!
!     *************** amplitude - phase notation ***********************
!
!             for i=1, ..., n
!
!             r(i) equals azero plus the sum from k=1 to k=kmax of
!
!                  alpha(k)*cos(k*(i-1)*2*pi/n+beta(k))
!
!             where
!
!                  alpha(k) = sqrt(a(k)*a(k)+b(k)*b(k))
!
!                  cos(beta(k))=a(k)/alpha(k)
!
!                  sin(beta(k))=-b(k)/alpha(k)
!
! **********************************************************************
!
!     subroutine sinti(n, wsave)
!
!     subroutine sinti initializes the array wsave which is used in
!     subroutine sint. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed.  the method
!             is most efficient when n+1 is a product of small primes.
!
!     output parameter
!
!     wsave   a work array with at least int(2.5*n+15) locations.
!             different wsave arrays are required for different values
!             of n. the contents of wsave must not be changed between
!             calls of sint.
!
! **********************************************************************
!
!     subroutine sint(n, x, wsave)
!
!     subroutine sint computes the discrete fourier sine transform
!     of an odd sequence x(i). the transform is defined below at
!     output parameter x.
!
!     sint is the unnormalized inverse of itself since a call of sint
!     followed by another call of sint will multiply the input sequence
!     x by 2*(n+1).
!
!     the array wsave which is used by subroutine sint must be
!     initialized by calling subroutine sinti(n, wsave).
!
!     input parameters
!
!     n       the length of the sequence to be transformed.  the method
!             is most efficient when n+1 is the product of small primes.
!
!     x       an array which contains the sequence to be transformed
!
!
!     wsave   a work array with dimension at least int(2.5*n+15)
!             in the program that calls sint. the wsave array must be
!             initialized by calling subroutine sinti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                  x(i)= the sum from k=1 to k=n
!
!                       2*x(k)*sin(k*i*pi/(n+1))
!
!                  a call of sint followed by another call of
!                  sint will multiply the sequence x by 2*(n+1).
!                  hence sint is the unnormalized inverse
!                  of itself.
!
!     wsave   contains initialization calculations which must not be
!             destroyed between calls of sint.
!
! **********************************************************************
!
!     subroutine costi(n, wsave)
!
!     subroutine costi initializes the array wsave which is used in
!     subroutine cost. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed.  the method
!             is most efficient when n-1 is a product of small primes.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             different wsave arrays are required for different values
!             of n. the contents of wsave must not be changed between
!             calls of cost.
!
! **********************************************************************
!
!     subroutine cost(n, x, wsave)
!
!     subroutine cost computes the discrete fourier cosine transform
!     of an even sequence x(i). the transform is defined below at output
!     parameter x.
!
!     cost is the unnormalized inverse of itself since a call of cost
!     followed by another call of cost will multiply the input sequence
!     x by 2*(n-1). the transform is defined below at output parameter x
!
!     the array wsave which is used by subroutine cost must be
!     initialized by calling subroutine costi(n, wsave).
!
!     input parameters
!
!     n       the length of the sequence x. n must be greater than 1.
!             the method is most efficient when n-1 is a product of
!             small primes.
!
!     x       an array which contains the sequence to be transformed
!
!     wsave   a work array which must be dimensioned at least 3*n+15
!             in the program that calls cost. the wsave array must be
!             initialized by calling subroutine costi(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                 x(i) = x(1)+(-1)**(i-1)*x(n)
!
!                  + the sum from k=2 to k=n-1
!
!                      2*x(k)*cos((k-1)*(i-1)*pi/(n-1))
!
!                  a call of cost followed by another call of
!                  cost will multiply the sequence x by 2*(n-1)
!                  hence cost is the unnormalized inverse
!                  of itself.
!
!     wsave   contains initialization calculations which must not be
!             destroyed between calls of cost.
!
! **********************************************************************
!
!     subroutine sinqi(n, wsave)
!
!     subroutine sinqi initializes the array wsave which is used in
!     both sinqf and sinqb. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed. the method
!             is most efficient when n is a product of small primes.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             the same work array can be used for both sinqf and sinqb
!             as long as n remains unchanged. different wsave arrays
!             are required for different values of n. the contents of
!             wsave must not be changed between calls of sinqf or sinqb.
!
! **********************************************************************
!
!     subroutine sinqf(n, x, wsave)
!
!     subroutine sinqf computes the fast fourier transform of quarter
!     wave data. that is , sinqf computes the coefficients in a sine
!     series representation with only odd wave numbers. the transform
!     is defined below at output parameter x.
!
!     sinqb is the unnormalized inverse of sinqf since a call of sinqf
!     followed by a call of sinqb will multiply the input sequence x
!     by 4*n.
!
!     the array wsave which is used by subroutine sinqf must be
!     initialized by calling subroutine sinqi(n, wsave).
!
!
!     input parameters
!
!     n       the length of the array x to be transformed.  the method
!             is most efficient when n is a product of small primes.
!
!     x       an array which contains the sequence to be transformed
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             in the program that calls sinqf. the wsave array must be
!             initialized by calling subroutine sinqi(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                  x(i) = (-1)**(i-1)*x(n)
!
!                     + the sum from k=1 to k=n-1 of
!
!                     2*x(k)*sin((2*i-1)*k*pi/(2*n))
!
!                  a call of sinqf followed by a call of
!                  sinqb will multiply the sequence x by 4*n.
!                  therefore sinqb is the unnormalized inverse
!                  of sinqf.
!
!     wsave   contains initialization calculations which must not
!             be destroyed between calls of sinqf or sinqb.
!
! **********************************************************************
!
!     subroutine sinqb(n, x, wsave)
!
!     subroutine sinqb computes the fast fourier transform of quarter
!     wave data. that is , sinqb computes a sequence from its
!     representation in terms of a sine series with odd wave numbers.
!     the transform is defined below at output parameter x.
!
!     sinqf is the unnormalized inverse of sinqb since a call of sinqb
!     followed by a call of sinqf will multiply the input sequence x
!     by 4*n.
!
!     the array wsave which is used by subroutine sinqb must be
!     initialized by calling subroutine sinqi(n, wsave).
!
!
!     input parameters
!
!     n       the length of the array x to be transformed.  the method
!             is most efficient when n is a product of small primes.
!
!     x       an array which contains the sequence to be transformed
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             in the program that calls sinqb. the wsave array must be
!             initialized by calling subroutine sinqi(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                  x(i)= the sum from k=1 to k=n of
!
!                    4*x(k)*sin((2k-1)*i*pi/(2*n))
!
!                  a call of sinqb followed by a call of
!                  sinqf will multiply the sequence x by 4*n.
!                  therefore sinqf is the unnormalized inverse
!                  of sinqb.
!
!     wsave   contains initialization calculations which must not
!             be destroyed between calls of sinqb or sinqf.
!
! **********************************************************************
!
!     subroutine cosqi(n, wsave)
!
!     subroutine cosqi initializes the array wsave which is used in
!     both cosqf and cosqb. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the array to be transformed.  the method
!             is most efficient when n is a product of small primes.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 3*n+15.
!             the same work array can be used for both cosqf and cosqb
!             as long as n remains unchanged. different wsave arrays
!             are required for different values of n. the contents of
!             wsave must not be changed between calls of cosqf or cosqb.
!
! **********************************************************************
!
!     subroutine cosqf(n, x, wsave)
!
!     subroutine cosqf computes the fast fourier transform of quarter
!     wave data. that is , cosqf computes the coefficients in a cosine
!     series representation with only odd wave numbers. the transform
!     is defined below at output parameter x
!
!     cosqf is the unnormalized inverse of cosqb since a call of cosqf
!     followed by a call of cosqb will multiply the input sequence x
!     by 4*n.
!
!     the array wsave which is used by subroutine cosqf must be
!     initialized by calling subroutine cosqi(n, wsave).
!
!
!     input parameters
!
!     n       the length of the array x to be transformed.  the method
!             is most efficient when n is a product of small primes.
!
!     x       an array which contains the sequence to be transformed
!
!     wsave   a work array which must be dimensioned at least 3*n+15
!             in the program that calls cosqf. the wsave array must be
!             initialized by calling subroutine cosqi(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                  x(i) = x(1) plus the sum from k=2 to k=n of
!
!                     2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))
!
!                  a call of cosqf followed by a call of
!                  cosqb will multiply the sequence x by 4*n.
!                  therefore cosqb is the unnormalized inverse
!                  of cosqf.
!
!     wsave   contains initialization calculations which must not
!             be destroyed between calls of cosqf or cosqb.
!
! **********************************************************************
!
!     subroutine cosqb(n, x, wsave)
!
!     subroutine cosqb computes the fast fourier transform of quarter
!     wave data. that is , cosqb computes a sequence from its
!     representation in terms of a cosine series with odd wave numbers.
!     the transform is defined below at output parameter x.
!
!     cosqb is the unnormalized inverse of cosqf since a call of cosqb
!     followed by a call of cosqf will multiply the input sequence x
!     by 4*n.
!
!     the array wsave which is used by subroutine cosqb must be
!     initialized by calling subroutine cosqi(n, wsave).
!
!
!     input parameters
!
!     n       the length of the array x to be transformed.  the method
!             is most efficient when n is a product of small primes.
!
!     x       an array which contains the sequence to be transformed
!
!     wsave   a work array that must be dimensioned at least 3*n+15
!             in the program that calls cosqb. the wsave array must be
!             initialized by calling subroutine cosqi(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!
!     output parameters
!
!     x       for i=1, ..., n
!
!                  x(i)= the sum from k=1 to k=n of
!
!                    4*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))
!
!                  a call of cosqb followed by a call of
!                  cosqf will multiply the sequence x by 4*n.
!                  therefore cosqf is the unnormalized inverse
!                  of cosqb.
!
!     wsave   contains initialization calculations which must not
!             be destroyed between calls of cosqb or cosqf.
!
! **********************************************************************
!
!     subroutine cffti(n, wsave)
!
!     subroutine cffti initializes the array wsave which is used in
!     both cfftf and cfftb. the prime factorization of n together with
!     a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 4*n+15
!             the same work array can be used for both cfftf and cfftb
!             as long as n remains unchanged. different wsave arrays
!             are required for different values of n. the contents of
!             wsave must not be changed between calls of cfftf or cfftb.
!
! **********************************************************************
!
!     subroutine cfftf(n, c, wsave)
!
!     subroutine cfftf computes the forward complex discrete fourier
!     transform (the fourier analysis). equivalently , cfftf computes
!     the fourier coefficients of a complex periodic sequence.
!     the transform is defined below at output parameter c.
!
!     the transform is not normalized. to obtain a normalized transform
!     the output must be divided by n. otherwise a call of cfftf
!     followed by a call of cfftb will multiply the sequence by n.
!
!     the array wsave which is used by subroutine cfftf must be
!     initialized by calling subroutine cffti(n, wsave).
!
!     input parameters
!
!
!     n      the length of the complex sequence c. the method is
!            more efficient when n is the product of small primes. n
!
!     c      a complex array of length n which contains the sequence
!
!     wsave   a real work array which must be dimensioned at least 4n+15
!             in the program that calls cfftf. the wsave array must be
!             initialized by calling subroutine cffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by cfftf and cfftb.
!
!     output parameters
!
!     c      for j=1, ..., n
!
!                c(j)=the sum from k=1, ..., n of
!
!                      c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)
!
!                            where i=sqrt(-1)
!
!     wsave   contains initialization calculations which must not be
!             destroyed between calls of subroutine cfftf or cfftb
!
! **********************************************************************
!
!     subroutine cfftb(n, c, wsave)
!
!     subroutine cfftb computes the backward complex discrete fourier
!     transform (the fourier synthesis). equivalently , cfftb computes
!     a complex periodic sequence from its fourier coefficients.
!     the transform is defined below at output parameter c.
!
!     a call of cfftf followed by a call of cfftb will multiply the
!     sequence by n.
!
!     the array wsave which is used by subroutine cfftb must be
!     initialized by calling subroutine cffti(n, wsave).
!
!     input parameters
!
!
!     n      the length of the complex sequence c. the method is
!            more efficient when n is the product of small primes.
!
!     c      a complex array of length n which contains the sequence
!
!     wsave   a real work array which must be dimensioned at least 4n+15
!             in the program that calls cfftb. the wsave array must be
!             initialized by calling subroutine cffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by cfftf and cfftb.
!
!     output parameters
!
!     c      for j=1, ..., n
!
!                c(j)=the sum from k=1, ..., n of
!
!                      c(k)*exp(i*(j-1)*(k-1)*2*pi/n)
!
!                            where i=sqrt(-1)
!
!     wsave   contains initialization calculations which must not be
!             destroyed between calls of subroutine cfftf or cfftb
!
module type_FastFourierTransform

    use spherepack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI, &
        TWO_PI, &
        HALF_PI

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: THREE = 3.0_wp
    real(wp), parameter :: FOUR = 4.0_wp
    real(wp), parameter :: FIVE = 5.0_wp
    real(wp), parameter :: SQRT_2 = sqrt(TWO)
    real(wp), parameter :: SQRT_3 = sqrt(THREE)
    real(wp), parameter :: SQRT_5 = sqrt(FIVE)

    type, public :: FastFourierTransform
    contains
        ! Type-bound procedures
        procedure, nopass, public :: rffti
        procedure, nopass, public :: rfftf
        procedure, nopass, public :: rfftb
        procedure, nopass, public :: ezffti
        procedure, nopass, public :: ezfftf
        procedure, nopass, public :: ezfftb
        procedure, nopass, public :: sinti
        procedure, nopass, public :: sint
        procedure, nopass, public :: costi
        procedure, nopass, public :: cost
        procedure, nopass, public :: sinqi
        procedure, nopass, public :: sinqf
        procedure, nopass, public :: sinqb
        procedure, nopass, public :: cosqi
        procedure, nopass, public :: cosqf
        procedure, nopass, public :: cosqb
        procedure, nopass, public :: cffti
        procedure, nopass, public :: cfftf
        procedure, nopass, public :: cfftb
    end type FastFourierTransform

contains

    subroutine ezfftf(n, r, azero, a, b, wsave)

        ! Dummy arguments
        integer(ip), intent(in)   :: n
        real(wp),    intent(out)  :: azero
        real(wp),    intent(in)   :: r(*)
        real(wp),    intent(out)  :: a(*)
        real(wp),    intent(out)  :: b(*)
        real(wp)                  :: wsave(*)

        ! Local variables
        integer(ip) :: ns2, ns2m
        real(wp)    :: cf, cfm

        if (3 > n) then
            if (n /= 2) then
                azero = r(1)
                return
            end if
            azero = HALF*(r(1)+r(2))
            a(1) = HALF*(r(1)-r(2))
            return
        end if

        wsave(:n) = r(:n)

        call rfftf(n, wsave, wsave(n+1))

        cf = TWO/n
        cfm = -cf
        azero = HALF*cf*wsave(1)
        ns2 =(n + 1)/2
        ns2m = ns2 - 1
        a(:ns2m) = cf*wsave(2:ns2m*2:2)
        b(:ns2m) = cfm*wsave(3:ns2m*2+1:2)

        if (mod(n, 2) /= 1) then
            a(ns2) = HALF*cf*wsave(n)
            b(ns2) = ZERO
        end if

    end subroutine ezfftf

    subroutine ezfftb(n, r, azero, a, b, wsave)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp), intent(in)    :: azero
        real(wp)                 :: r(*)
        real(wp), intent(in)    :: a(*)
        real(wp), intent(in)    :: b(*)
        real(wp)                 :: wsave(*)

        ! Local variables

        integer(ip) :: ns2


        if (3 > n) then
            if (n /= 2) then
                r(1) = azero
                return
            end if
            r(1) = azero + a(1)
            r(2) = azero - a(1)
            return
        end if

        ns2 =(n - 1)/2
        r(2:ns2*2:2) = HALF*a(:ns2)
        r(3:ns2*2+1:2) = -HALF*b(:ns2)
        r(1) = azero

        if (mod(n, 2) == 0) r(n) = a(ns2+1)

        call rfftb(n, r, wsave(n+1))

    end subroutine ezfftb

    subroutine ezffti(n, wsave)

        ! Dummy arguments
        integer(ip) :: n
        real(wp)    :: wsave(*)

        if (n > 1) then
            call real_periodic_initialization_lower_utility_routine(n, wsave(2*n+1), wsave(3*n+1))
        end if

    end subroutine ezffti

    subroutine real_periodic_initialization_lower_utility_routine(n, wa, ifac)

        ! Dummy arguments
        integer(ip),  intent(in)     :: n
        real(wp),     intent(inout)  :: ifac(*)
        real(wp),     intent(inout)  :: wa(*)

        ! Local variables
        integer(ip), parameter :: NTRYH(*) = [4, 2, 3, 5]
        integer(ip) :: nl, nf, j, ntry, nq, nr, i, is, nfm1
        integer(ip) :: l1, k1, iip, l2, ido, iipm, ii
        real(wp) :: argh, arg1, ch1, sh1, dch1, dsh1, temp

        ! Initialize
        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorization_loop: do

            j = j + 1

            if (j <= 4) then
                ntry = NTRYH(j)
            else
                ntry = ntry + 2
            end if

            inner_loop: do

                nq = nl/ntry
                nr = nl - ntry*nq

                if (nr /= 0) cycle factorization_loop

                nf = nf + 1
                ifac(nf+2) = ntry
                nl = nq

                if (ntry == 2) then
                    if (nf /= 1) then
                        ifac(nf+2:4:(-1)) = ifac(nf+1:3:(-1))
                        ifac(3) = 2
                    end if
                end if

                if (nl /= 1) cycle inner_loop

                exit inner_loop
            end do inner_loop
            exit factorization_loop
        end do factorization_loop

        ifac(1) = n
        ifac(2) = nf
        argh = TWO_PI/n
        is = 0
        nfm1 = nf - 1
        l1 = 1

        if (nfm1 /= 0) then
            do k1 = 1, nfm1
                iip = int(ifac(k1+2), kind=ip)
                l2 = l1*iip
                ido = n/l2
                iipm = iip - 1
                arg1 = real(l1, kind=wp)*argh
                ch1 = ONE
                sh1 = ZERO
                dch1 = cos(arg1)
                dsh1 = sin(arg1)
                do j = 1, iipm
                    temp = dch1*ch1 - dsh1*sh1
                    sh1 = dch1*sh1 + dsh1*ch1
                    ch1 = temp
                    i = is + 2
                    wa(i-1) = ch1
                    wa(i) = sh1
                    if (ido >= 5) then
                        do ii = 5, ido, 2
                            i = i + 2
                            wa(i-1) = ch1*wa(i-3) - sh1*wa(i-2)
                            wa(i) = ch1*wa(i-2) + sh1*wa(i-3)
                        end do
                    end if
                    is = is + ido
                end do
                l1 = l2
            end do
        end if

    end subroutine real_periodic_initialization_lower_utility_routine

    subroutine costi(n, wsave)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp)                 :: wsave(*)

        ! Local variables

        integer(ip)  :: nm1, np1, ns2, k, kc
        real(wp)     :: dt, fk


        if (n > 3) then
            nm1 = n - 1
            np1 = n + 1
            ns2 = n/2
            dt = pi/nm1
            fk = ZERO

            do k = 2, ns2
                kc = np1 - k
                fk = fk + ONE
                wsave(k) = TWO * sin(fk*dt)
                wsave(kc) = TWO * cos(fk*dt)
            end do

            call rffti(nm1, wsave(n+1))

        end if

    end subroutine costi

    subroutine cost(n, x, wsave)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp) :: x(*)
        real(wp) :: wsave(*)

        ! Local variables

        integer(ip) :: nm1, np1, ns2, k, kc, modn, i
        real(wp) :: x1h, x1p3, tx2, c1, t1, t2, xim2, xi

        !
        nm1 = n - 1
        np1 = n + 1
        ns2 = n/2

        if (3 < n) then
            if (3 > n) then
                x1h = x(1) + x(2)
                x(2) = x(1) - x(2)
                x(1) = x1h
                return
            end if
            if (n <= 3) then
                x1p3 = x(1) + x(3)
                tx2 = x(2) + x(2)
                x(2) = x(1) - x(3)
                x(1) = x1p3 + tx2
                x(3) = x1p3 - tx2
                return
            end if
            c1 = x(1) - x(n)
            x(1) = x(1) + x(n)
            do k = 2, ns2
                kc = np1 - k
                t1 = x(k) + x(kc)
                t2 = x(k) - x(kc)
                c1 = c1 + wsave(kc)*t2
                t2 = wsave(k)*t2
                x(k) = t1 - t2
                x(kc) = t1 + t2
            end do

            modn = mod(n, 2)

            if (modn /= 0) x(ns2+1) = x(ns2+1) + x(ns2+1)

            call rfftf(nm1, x, wsave(n+1))

            xim2 = x(2)
            x(2) = c1

            do i = 4, n, 2
                xi = x(i)
                x(i) = x(i-2) - x(i-1)
                x(i-1) = xim2
                xim2 = xi
            end do

            if (modn /= 0) x(n) = xim2

        end if

    end subroutine cost

    subroutine sinti(n, wsave)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp)                 :: wsave(*)

        ! Local variables

        integer(ip) :: ns2, np1, k
        real(wp)    :: dt


        if (n > 1) then
            ns2 = n/2
            np1 = n + 1
            dt = PI/np1

            do k = 1, ns2
                wsave(k) = TWO * sin(k*dt)
            end do

            call rffti(np1, wsave(ns2+1))
        end if

    end subroutine sinti

    subroutine sint(n, x, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: x(*)
        real(wp) :: wsave(*)

        ! Local variables

        integer(ip) :: np1, iw1, iw2, iw3


        np1 = n + 1
        iw1 = n/2 + 1
        iw2 = iw1 + np1
        iw3 = iw2 + np1

        call sint_lower_utility_routine(n, x, wsave, wsave(iw1), wsave(iw2), wsave(iw3))

    end subroutine sint

    subroutine sint_lower_utility_routine(n, war, was, xh, x, ifac)

        ! Dummy arguments

        integer(ip), intent(in)  :: n
        real(wp)              :: ifac(*)
        real(wp)              :: war(*)
        real(wp), intent(in) :: was(*)
        real(wp)              :: xh(*)
        real(wp)              :: x(*)

        ! Local variables

        integer(ip)         :: i, np1, ns2, k, kc, modn
        real(wp)            :: temp, t1, t2


        xh(:n) = war(:n)
        war(:n) = x(:n)

        if (3 > n) then
            if (n /= 2) then
                xh(1) = xh(1) + xh(1)
                x(:n) = war(:n)
                war(:n) = xh(:n)
                return
            end if
            temp = SQRT_3*(xh(1)+xh(2))
            xh(2) = SQRT_3*(xh(1)-xh(2))
            xh(1) = temp
            x(:n) = war(:n)
            war(:n) = xh(:n)
            return
        end if

        np1 = n + 1
        ns2 = n/2
        x(1) = ZERO

        do k = 1, ns2
            kc = np1 - k
            t1 = xh(k) - xh(kc)
            t2 = was(k)*(xh(k)+xh(kc))
            x(k+1) = t1 + t2
            x(kc+1) = t2 - t1
        end do

        modn = mod(n, 2)

        if (modn /= 0) x(ns2+2) = FOUR * xh(ns2+1)

        call real_forward_lower_utility_routine(np1, x, xh, war, ifac)

        xh(1) = HALF*x(1)

        do i = 3, n, 2
            xh(i-1) = -x(i)
            xh(i) = xh(i-2) + x(i-1)
        end do

        if (modn == 0) xh(n) = -x(n+1)

        x(:n) = war(:n)
        war(:n) = xh(:n)

    end subroutine sint_lower_utility_routine

    subroutine cosqi(n, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: wsave(*)

        ! Local variables

        integer(ip) :: k ! Counter
        real(wp)    :: fk, dt


        dt = HALF_PI/n

        fk = ZERO
        do k = 1, n
            fk = fk + ONE
            wsave(k) = cos(fk*dt)
        end do

        call rffti(n, wsave(n+1))

    end subroutine cosqi

    subroutine cosqf(n, x, wsave)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp)                 :: x(*)
        real(wp)                 :: wsave(*)

        ! Local variables

        real(wp) :: temp


        if (3 < n) then
            if (n > 2) call cosqf_lower_utility_routine(n, x, wsave, wsave(n+1))
            temp = SQRT_2*x(2)
            x(2) = x(1) - temp
            x(1) = x(1) + temp
        end if

    end subroutine cosqf

    subroutine cosqf_lower_utility_routine(n, x, w, xh)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp)                 :: x(*)
        real(wp),    intent(in) :: w(*)
        real(wp)                 :: xh(*)

        ! Local variables

        integer(ip) :: ns2, np2, k, kc, modn, i
        real(wp)    :: temp


        ns2 =(n + 1)/2
        np2 = n + 2

        do k = 2, ns2
            kc = np2 - k
            xh(k) = x(k) + x(kc)
            xh(kc) = x(k) - x(kc)
        end do

        modn = mod(n, 2)

        if (modn == 0) xh(ns2+1) = x(ns2+1) + x(ns2+1)

        do k = 2, ns2
            kc = np2 - k
            x(k) = w(k-1)*xh(kc) + w(kc-1)*xh(k)
            x(kc) = w(k-1)*xh(k) - w(kc-1)*xh(kc)
        end do

        if (modn == 0) x(ns2+1) = w(ns2)*xh(ns2+1)

        call rfftf(n, x, xh)

        do i = 3, n, 2
            temp = x(i-1) - x(i)
            x(i) = x(i-1) + x(i)
            x(i-1) = temp
        end do

    end subroutine cosqf_lower_utility_routine




    subroutine cosqb(n, x, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: x(*)
        real(wp) :: wsave(*)

        ! Local variables

        real(wp), parameter :: TWO_SQRT_2 = TWO * sqrt(TWO) ! 2.82842712474619
        real(wp)            :: temp


        if (3 > n) then
            if (n /= 2) then
                x(1) = FOUR * x(1)
            else
                temp = FOUR * (x(1)+x(2))
                x(2) = TWO_SQRT_2*(x(1)-x(2))
                x(1) = temp
            end if
        end if

        call cosqb_lower_utility_routine(n, x, wsave, wsave(n+1))

    end subroutine cosqb


    subroutine cosqb_lower_utility_routine(n, x, w, xh)

        ! Dummy arguments

        integer(ip)               :: n
        real(wp)              :: x(*)
        real(wp), intent(in) :: w(*)
        real(wp)              :: xh(*)

        ! Local variables

        integer(ip) :: ns2, np2, i, modn, k, kc
        real(wp)    :: xim1


        ns2 =(n + 1)/2
        np2 = n + 2

        do i = 3, n, 2
            xim1 = x(i-1) + x(i)
            x(i) = x(i) - x(i-1)
            x(i-1) = xim1
        end do

        x(1) = TWO*x(1)

        modn = mod(n, 2)

        if (modn == 0) x(n) = x(n) + x(n)

        call rfftb(n, x, xh)

        do k = 2, ns2
            kc = np2 - k
            xh(k) = w(k-1)*x(kc) + w(kc-1)*x(k)
            xh(kc) = w(k-1)*x(k) - w(kc-1)*x(kc)
        end do

        if (modn == 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))

        do k = 2, ns2
            kc = np2 - k
            x(k) = xh(k) + xh(kc)
            x(kc) = xh(k) - xh(kc)
        end do

        x(1) = TWO*x(1)

    end subroutine cosqb_lower_utility_routine

    subroutine sinqi(n, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: wsave(*)


        call cosqi(n, wsave)

    end subroutine sinqi

    subroutine sinqf(n, x, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: x(*)
        real(wp) :: wsave(*)

        ! Local variables

        integer(ip) :: ns2, k, kc
        real(wp) :: xhold


        if (n /= 1) then

            ns2 = n/2

            do k = 1, ns2
                kc = n - k
                xhold = x(k)
                x(k) = x(kc+1)
                x(kc+1) = xhold
            end do

            call cosqf(n, x, wsave)

            x(2:n:2) = -x(2:n:2)

        end if

    end subroutine sinqf

    subroutine sinqb(n, x, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: x(*)
        real(wp) :: wsave(*)

        ! Local variables

        integer(ip) :: ns2, k, kc
        real(wp) :: xhold


        if (n <= 1) then
            x(1) = FOUR * x(1)
        else
            ns2 = n/2
            x(2:n:2) = -x(2:n:2)

            call cosqb(n, x, wsave)

            do k = 1, ns2
                kc = n - k
                xhold = x(k)
                x(k) = x(kc+1)
                x(kc+1) = xhold
            end do
        end if

    end subroutine sinqb

    subroutine cffti(n, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: wsave(*)

        ! Dummy arguments

        integer(ip) :: iw1, iw2


        if (n /= 1) then
            iw1 = 2*n + 1
            iw2 = iw1 + 2*n
            call complex_initialization_lower_utility_routine(n, wsave(iw1), wsave(iw2))
        end if

    end subroutine cffti

    subroutine complex_initialization_lower_utility_routine(n, wa, ifac)

        ! Dummy arguments

        integer(ip), intent(in)      :: n
        real(wp), intent(inout)  :: ifac(*)
        real(wp), intent(inout)  :: wa(*)

        ! Local variables

        integer(ip), parameter  :: NTRYH(*) = [3, 4, 2, 5]
        integer(ip)             :: nl, nf, j, ntry, nq, nr
        integer(ip)             :: i, l1, k1, iip, ld, l2, ido
        integer(ip)             :: idot, iipm, i1, ii
        real(wp)                :: argh, fi, argld, arg


        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorization_loop: do

            j = j + 1

            if (j - 4 <= 0) then
                ntry = NTRYH(j)
            else
                ntry = ntry + 2
            end if

            inner_loop: do

                nq = nl/ntry
                nr = nl - ntry*nq

                if (nr /= 0) cycle factorization_loop

                nf = nf + 1
                ifac(nf+2) = ntry
                nl = nq

                if (ntry == 2) then
                    if (nf /= 1) then
                        ifac(nf+2:4:(-1)) = ifac(nf+1:3:(-1))
                        ifac(3) = 2
                    end if
                end if

                if (nl /= 1) cycle inner_loop
                exit inner_loop
            end do inner_loop
            exit factorization_loop
        end do factorization_loop

        ifac(1) = n
        ifac(2) = nf
        argh = TWO_PI/n
        i = 2
        l1 = 1
        do k1 = 1, nf
            iip = int(ifac(k1+2), kind=ip)
            ld = 0
            l2 = l1*iip
            ido = n/l2
            idot = ido + ido + 2
            iipm = iip - 1
            set_workspace_loop: do j = 1, iipm
                i1 = i
                wa(i-1) = ONE
                wa(i) = ZERO
                ld = ld + l1
                fi = ZERO
                argld = real(ld, kind=wp)*argh
                do ii = 4, idot, 2
                    i = i + 2
                    fi = fi + ONE
                    arg = fi*argld
                    wa(i-1) = cos(arg)
                    wa(i) = sin(arg)
                end do

                if (iip <= 5) cycle set_workspace_loop

                wa(i1-1) = wa(i-1)
                wa(i1) = wa(i)

            end do set_workspace_loop
            l1 = l2
        end do

    end subroutine complex_initialization_lower_utility_routine

    subroutine cfftb(n, c, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: c(*)
        real(wp) :: wsave(*)

        ! Dummy arguments

        integer(ip) :: iw1, iw2


        if (n /= 1) then
            iw1 = 2*n + 1
            iw2 = iw1 + 2*n
            call complex_backward_lower_utility_routine(n, c, wsave, wsave(iw1), wsave(iw2))
        end if

    end subroutine cfftb

    subroutine complex_backward_lower_utility_routine(n, c, ch, wa, ifac)

        ! Dummy arguments

        integer(ip), intent(in) :: n
        real(wp), intent(in) :: ifac(*)
        real(wp) :: c(*)
        real(wp) :: ch(*)
        real(wp) :: wa(*)

        ! Local variables

        integer(ip) :: nf, na, l1, iw, k1, iip, l2, ido
        integer(ip) :: idot, idl1, ix2, ix3, ix4, nac, n2


        nf = int(ifac(2), kind=ip)
        na = 0
        l1 = 1
        iw = 1
        do k1 = 1, nf
            iip = int(ifac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idot = ido + ido
            idl1 = idot*l1

            select case (iip)
                case (2)
                    if (na == 0) then
                        call complex_backward_pass_2(idot, l1, c, ch, wa(iw))
                    else
                        call complex_backward_pass_2(idot, l1, ch, c, wa(iw))
                    end if
                    na = 1 - na
                case (3)
                    ix2 = iw + idot
                    if (na == 0) then
                        call complex_backward_pass_3(idot, l1, c, ch, wa(iw), wa(ix2))
                    else
                        call complex_backward_pass_3(idot, l1, ch, c, wa(iw), wa(ix2))
                    end if
                    na = 1 - na
                case (4)
                    ix2 = iw + idot
                    ix3 = ix2 + idot
                    if (na == 0) then
                        call complex_backward_pass_4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
                    else
                        call complex_backward_pass_4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
                    end if
                    na = 1 - na
                case (5)
                    ix2 = iw + idot
                    ix3 = ix2 + idot
                    ix4 = ix3 + idot
                    if (na == 0) then
                        call complex_backward_pass_5(idot, l1, c, ch, wa(iw), wa(ix2), &
                            wa(ix3), wa(ix4))
                    else
                        call complex_backward_pass_5(idot, l1, ch, c, wa(iw), wa(ix2), &
                            wa(ix3), wa(ix4))
                    end if
                    na = 1 - na
                case default
                    if (na == 0) then
                        call complex_backward_pass_n(nac, idot, iip, l1, idl1, c, c, c, ch, &
                            ch, wa(iw))
                    else
                        call complex_backward_pass_n(nac, idot, iip, l1, idl1, ch, ch, ch, &
                            c, c, wa(iw))
                    end if
                    if (nac /= 0) then
                        na = 1 - na
                    end if
            end select
            l1 = l2
            iw = iw + (iip - 1)*idot
        end do

        if (na /= 0) then
            n2 = 2*n
            c(:n2) = ch(:n2)
        end if

    end subroutine complex_backward_lower_utility_routine

    pure subroutine complex_backward_pass_2(ido, l1, cc, ch, wa1)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 2, l1)
        real(wp), intent(out) :: ch(ido, l1, 2)
        real(wp), intent(in) :: wa1(1)

        ! Local variables

        integer(ip) :: k, i
        real(wp) :: tr2, ti2


        if(ido <= 2) then
            ch(1, :, 1) = cc(1, 1, :) + cc(1, 2, :)
            ch(1, :, 2) = cc(1, 1, :) - cc(1, 2, :)
            ch(2, :, 1) = cc(2, 1, :) + cc(2, 2, :)
            ch(2, :, 2) = cc(2, 1, :) - cc(2, 2, :)
        else
            do k = 1, l1
                do i = 2, ido, 2
                    ch(i-1, k, 1) = cc(i-1, 1, k) + cc(i-1, 2, k)
                    tr2 = cc(i-1, 1, k) - cc(i-1, 2, k)
                    ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
                    ti2 = cc(i, 1, k) - cc(i, 2, k)
                    ch(i, k, 2) = wa1(i-1)*ti2 + wa1(i)*tr2
                    ch(i-1, k, 2) = wa1(i-1)*tr2 - wa1(i)*ti2
                end do
            end do
        end if

    end subroutine complex_backward_pass_2

    pure subroutine complex_backward_pass_3(ido, l1, cc, ch, wa1, wa2)

        ! Dummy arguments

        integer(ip), intent(in)  :: ido
        integer(ip), intent(in)  :: l1
        real(wp),    intent(in)  :: cc(ido, 3, l1)
        real(wp),    intent(out) :: ch(ido, l1, 3)
        real(wp),    intent(in)  :: wa1(*)
        real(wp),    intent(in)  :: wa2(*)

        ! Local variables

        integer(ip)         :: k, i
        real(wp), parameter :: taur = -HALF
        real(wp), parameter :: taui = SQRT_3/2 ! 0.866025403784439
        real(wp)            :: tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3


        select case(ido)
            case (2)
                do k = 1, l1
                    tr2 = cc(1, 2, k) + cc(1, 3, k)
                    cr2 = cc(1, 1, k) + taur*tr2
                    ch(1, k, 1) = cc(1, 1, k) + tr2
                    ti2 = cc(2, 2, k) + cc(2, 3, k)
                    ci2 = cc(2, 1, k) + taur*ti2
                    ch(2, k, 1) = cc(2, 1, k) + ti2
                    cr3 = taui*(cc(1, 2, k)-cc(1, 3, k))
                    ci3 = taui*(cc(2, 2, k)-cc(2, 3, k))
                    ch(1, k, 2) = cr2 - ci3
                    ch(1, k, 3) = cr2 + ci3
                    ch(2, k, 2) = ci2 + cr3
                    ch(2, k, 3) = ci2 - cr3
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        tr2 = cc(i-1, 2, k) + cc(i-1, 3, k)
                        cr2 = cc(i-1, 1, k) + taur*tr2
                        ch(i-1, k, 1) = cc(i-1, 1, k) + tr2
                        ti2 = cc(i, 2, k) + cc(i, 3, k)
                        ci2 = cc(i, 1, k) + taur*ti2
                        ch(i, k, 1) = cc(i, 1, k) + ti2
                        cr3 = taui*(cc(i-1, 2, k)-cc(i-1, 3, k))
                        ci3 = taui*(cc(i, 2, k)-cc(i, 3, k))
                        dr2 = cr2 - ci3
                        dr3 = cr2 + ci3
                        di2 = ci2 + cr3
                        di3 = ci2 - cr3
                        ch(i, k, 2) = wa1(i-1)*di2 + wa1(i)*dr2
                        ch(i-1, k, 2) = wa1(i-1)*dr2 - wa1(i)*di2
                        ch(i, k, 3) = wa2(i-1)*di3 + wa2(i)*dr3
                        ch(i-1, k, 3) = wa2(i-1)*dr3 - wa2(i)*di3
                    end do
                end do
        end select

    end subroutine complex_backward_pass_3

    pure subroutine complex_backward_pass_4(ido, l1, cc, ch, wa1, wa2, wa3)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 4, l1)
        real(wp), intent(out) :: ch(ido, l1, 4)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp) :: ti1, ti2, tr4, ti3, tr1, tr2, ti4
        real(wp) :: tr3, cr3, ci3, cr2, cr4, ci2, ci4


        select case(ido)
            case (2)
                do k = 1, l1
                    ti1 = cc(2, 1, k) - cc(2, 3, k)
                    ti2 = cc(2, 1, k) + cc(2, 3, k)
                    tr4 = cc(2, 4, k) - cc(2, 2, k)
                    ti3 = cc(2, 2, k) + cc(2, 4, k)
                    tr1 = cc(1, 1, k) - cc(1, 3, k)
                    tr2 = cc(1, 1, k) + cc(1, 3, k)
                    ti4 = cc(1, 2, k) - cc(1, 4, k)
                    tr3 = cc(1, 2, k) + cc(1, 4, k)
                    ch(1, k, 1) = tr2 + tr3
                    ch(1, k, 3) = tr2 - tr3
                    ch(2, k, 1) = ti2 + ti3
                    ch(2, k, 3) = ti2 - ti3
                    ch(1, k, 2) = tr1 + tr4
                    ch(1, k, 4) = tr1 - tr4
                    ch(2, k, 2) = ti1 + ti4
                    ch(2, k, 4) = ti1 - ti4
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        ti1 = cc(i, 1, k) - cc(i, 3, k)
                        ti2 = cc(i, 1, k) + cc(i, 3, k)
                        ti3 = cc(i, 2, k) + cc(i, 4, k)
                        tr4 = cc(i, 4, k) - cc(i, 2, k)
                        tr1 = cc(i-1, 1, k) - cc(i-1, 3, k)
                        tr2 = cc(i-1, 1, k) + cc(i-1, 3, k)
                        ti4 = cc(i-1, 2, k) - cc(i-1, 4, k)
                        tr3 = cc(i-1, 2, k) + cc(i-1, 4, k)
                        ch(i-1, k, 1) = tr2 + tr3
                        cr3 = tr2 - tr3
                        ch(i, k, 1) = ti2 + ti3
                        ci3 = ti2 - ti3
                        cr2 = tr1 + tr4
                        cr4 = tr1 - tr4
                        ci2 = ti1 + ti4
                        ci4 = ti1 - ti4
                        ch(i-1, k, 2) = wa1(i-1)*cr2 - wa1(i)*ci2
                        ch(i, k, 2) = wa1(i-1)*ci2 + wa1(i)*cr2
                        ch(i-1, k, 3) = wa2(i-1)*cr3 - wa2(i)*ci3
                        ch(i, k, 3) = wa2(i-1)*ci3 + wa2(i)*cr3
                        ch(i-1, k, 4) = wa3(i-1)*cr4 - wa3(i)*ci4
                        ch(i, k, 4) = wa3(i-1)*ci4 + wa3(i)*cr4
                    end do
                end do
        end select

    end subroutine complex_backward_pass_4

    pure subroutine complex_backward_pass_5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 5, l1)
        real(wp), intent(out) :: ch(ido, l1, 5)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)
        real(wp), intent(in) :: wa4(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp), parameter :: tr11 = (SQRT_5 - ONE)/4 ! 0.309016994374947
        real(wp), parameter :: ti11 = (sqrt(ONE/(FIVE + SQRT_5)))/2 ! 0.951056516295154
        real(wp), parameter :: tr12 = (-ONE - SQRT_5)/4 ! -.809016994374947
        real(wp), parameter :: ti12 = sqrt( FIVE/(TWO*(FIVE + SQRT_5))) ! 0.587785252292473
        real(wp) ::  ti5, ti2, ti4, ti3, tr5, tr2, tr4
        real(wp) :: tr3, cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4, dr3, dr4, di3
        real(wp) :: di4, dr5, dr2, di5, di2


        select case(ido)
            case (2)
                do k = 1, l1
                    ti5 = cc(2, 2, k) - cc(2, 5, k)
                    ti2 = cc(2, 2, k) + cc(2, 5, k)
                    ti4 = cc(2, 3, k) - cc(2, 4, k)
                    ti3 = cc(2, 3, k) + cc(2, 4, k)
                    tr5 = cc(1, 2, k) - cc(1, 5, k)
                    tr2 = cc(1, 2, k) + cc(1, 5, k)
                    tr4 = cc(1, 3, k) - cc(1, 4, k)
                    tr3 = cc(1, 3, k) + cc(1, 4, k)
                    ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
                    ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
                    cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
                    ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
                    cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
                    ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
                    cr5 = ti11*tr5 + ti12*tr4
                    ci5 = ti11*ti5 + ti12*ti4
                    cr4 = ti12*tr5 - ti11*tr4
                    ci4 = ti12*ti5 - ti11*ti4
                    ch(1, k, 2) = cr2 - ci5
                    ch(1, k, 5) = cr2 + ci5
                    ch(2, k, 2) = ci2 + cr5
                    ch(2, k, 3) = ci3 + cr4
                    ch(1, k, 3) = cr3 - ci4
                    ch(1, k, 4) = cr3 + ci4
                    ch(2, k, 4) = ci3 - cr4
                    ch(2, k, 5) = ci2 - cr5
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        ti5 = cc(i, 2, k) - cc(i, 5, k)
                        ti2 = cc(i, 2, k) + cc(i, 5, k)
                        ti4 = cc(i, 3, k) - cc(i, 4, k)
                        ti3 = cc(i, 3, k) + cc(i, 4, k)
                        tr5 = cc(i-1, 2, k) - cc(i-1, 5, k)
                        tr2 = cc(i-1, 2, k) + cc(i-1, 5, k)
                        tr4 = cc(i-1, 3, k) - cc(i-1, 4, k)
                        tr3 = cc(i-1, 3, k) + cc(i-1, 4, k)
                        ch(i-1, k, 1) = cc(i-1, 1, k) + tr2 + tr3
                        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
                        cr2 = cc(i-1, 1, k) + tr11*tr2 + tr12*tr3
                        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
                        cr3 = cc(i-1, 1, k) + tr12*tr2 + tr11*tr3
                        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
                        cr5 = ti11*tr5 + ti12*tr4
                        ci5 = ti11*ti5 + ti12*ti4
                        cr4 = ti12*tr5 - ti11*tr4
                        ci4 = ti12*ti5 - ti11*ti4
                        dr3 = cr3 - ci4
                        dr4 = cr3 + ci4
                        di3 = ci3 + cr4
                        di4 = ci3 - cr4
                        dr5 = cr2 + ci5
                        dr2 = cr2 - ci5
                        di5 = ci2 - cr5
                        di2 = ci2 + cr5
                        ch(i-1, k, 2) = wa1(i-1)*dr2 - wa1(i)*di2
                        ch(i, k, 2) = wa1(i-1)*di2 + wa1(i)*dr2
                        ch(i-1, k, 3) = wa2(i-1)*dr3 - wa2(i)*di3
                        ch(i, k, 3) = wa2(i-1)*di3 + wa2(i)*dr3
                        ch(i-1, k, 4) = wa3(i-1)*dr4 - wa3(i)*di4
                        ch(i, k, 4) = wa3(i-1)*di4 + wa3(i)*dr4
                        ch(i-1, k, 5) = wa4(i-1)*dr5 - wa4(i)*di5
                        ch(i, k, 5) = wa4(i-1)*di5 + wa4(i)*dr5
                    end do
                end do
        end select

    end subroutine complex_backward_pass_5

    subroutine complex_backward_pass_n(nac, ido, iip, l1, idl1, cc, c1, c2, ch, ch2, wa)

        ! Dummy arguments

        integer(ip), intent(out) :: nac
        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: iip
        integer(ip), intent(in) :: l1
        integer(ip), intent(in) :: idl1
        real(wp), intent(in) :: cc(ido, iip, l1)
        real(wp), intent(out) :: c1(ido, l1, iip)
        real(wp), intent(inout)  :: c2(idl1, iip)
        real(wp), intent(inout)  :: ch(ido, l1, iip)
        real(wp), intent(inout)  :: ch2(idl1, iip)
        real(wp), intent(in) :: wa(*)

        ! Local variables

        integer(ip) :: idot, nt, iipp2, iipph, idp, j, jc, k
        integer(ip) :: i, idl, inc, l, lc, idlj, idij, idj
        real(wp) :: war, wai


        idot = ido/2
        nt = iip*idl1
        iipp2 = iip + 2
        iipph = (iip + 1)/2
        idp = iip*ido

        if (l1 <= ido) then
            do j = 2, iipph
                jc = iipp2 - j
                ch(:, :, j) = cc(:, j, :) + cc(:, jc, :)
                ch(:, :, jc) = cc(:, j, :) - cc(:, jc, :)
            end do
            ch(:, :, 1) = cc(:, 1, :)
        else
            do j = 2, iipph
                jc = iipp2 - j
                ch(:, :, j) = cc(:, j, :) + cc(:, jc, :)
                ch(:, :, jc) = cc(:, j, :) - cc(:, jc, :)
            end do
            ch(:, :, 1) = cc(:, 1, :)
        end if

        idl = 2 - ido
        inc = 0

        do l = 2, iipph
            lc = iipp2 - l
            idl = idl + ido
            c2(:, l) = ch2(:, 1) + wa(idl-1)*ch2(:, 2)
            c2(:, lc) = wa(idl)*ch2(:, iip)
            idlj = idl
            inc = inc + ido
            do j = 3, iipph
                jc = iipp2 - j
                idlj = idlj + inc
                if (idlj > idp) idlj = idlj - idp
                war = wa(idlj-1)
                wai = wa(idlj)
                c2(:, l) = c2(:, l) + war*ch2(:, j)
                c2(:, lc) = c2(:, lc) + wai*ch2(:, jc)
            end do
        end do

        do j = 2, iipph
            ch2(:, 1) = ch2(:, 1) + ch2(:, j)
        end do

        do j = 2, iipph
            jc = iipp2 - j
            ch2(:idl1-1:2, j) = c2(:idl1-1:2, j) - c2(2:idl1:2, jc)
            ch2(:idl1-1:2, jc) = c2(:idl1-1:2, j) + c2(2:idl1:2, jc)
            ch2(2:idl1:2, j) = c2(2:idl1:2, j) + c2(:idl1-1:2, jc)
            ch2(2:idl1:2, jc) = c2(2:idl1:2, j) - c2(:idl1-1:2, jc)
        end do

        nac = 1

        if(ido /= 2) then
            nac = 0
            c2(:, 1) = ch2(:, 1)
            c1(1, :, 2:iip) = ch(1, :, 2:iip)
            c1(2, :, 2:iip) = ch(2, :, 2:iip)

            if(idot <= l1) then
                idij = 0
                do j = 2, iip
                    idij = idij + 2
                    do i = 4, ido, 2
                        idij = idij + 2
                        c1(i-1, :, j) = wa(idij-1)*ch(i-1, :, j) - wa(idij)*ch(i, :, j)
                        c1(i, :, j) = wa(idij-1)*ch(i, :, j) + wa(idij)*ch(i-1, :, j)
                    end do
                end do
            else
                idj = 2 - ido
                do j = 2, iip
                    idj = idj + ido
                    do k = 1, l1
                        idij = idj
                        c1(3:ido-1:2, k, j) = wa(idij+1:ido-3+idij:2)*ch(3:ido-1:2, k, j &
                            ) - wa(idij+2:ido-2+idij:2)*ch(4:ido:2, k, j)
                        c1(4:ido:2, k, j) = wa(idij+1:ido-3+idij:2)*ch(4:ido:2, k, j) + &
                            wa(idij+2:ido-2+idij:2)*ch(3:ido-1:2, k, j)
                    end do
                end do
            end if
        end if

    end subroutine complex_backward_pass_n

    subroutine cfftf(n, c, wsave)

        ! Dummy arguments

        integer(ip), intent(in)      :: n
        real(wp),    intent(inout)   :: c(*)
        real(wp),    intent(inout)   :: wsave(*)

        ! Local variables

        integer(ip) :: iw1, iw2


        if (n /= 1) then
            iw1 = 2*n + 1
            iw2 = iw1 + 2*n
            call complex_forward_lower_utility_routine(n, c, wsave, wsave(iw1), wsave(iw2))
        end if

    end subroutine cfftf

    subroutine complex_forward_lower_utility_routine(n, c, ch, wa, ifac)

        ! Dummy arguments

        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: c(*)
        real(wp),    intent(inout)  :: ch(*)
        real(wp),    intent(in)     :: wa(*)
        real(wp),    intent(in)     :: ifac(*)

        ! Local variables

        integer(ip) :: nf, na, l1, iw, k1, iip, l2, ido
        integer(ip) :: idot, idl1, ix2, ix3, ix4, nac, n2


        nf = int(ifac(2), kind=ip)
        na = 0
        l1 = 1
        iw = 1

        do k1 = 1, nf
            iip = int(ifac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idot = ido + ido
            idl1 = idot*l1

            select case (iip)
                case (2)
                    if (na == 0) then
                        call complex_forward_pass_2(idot, l1, c, ch, wa(iw))
                    else
                        call complex_forward_pass_2(idot, l1, ch, c, wa(iw))
                    end if
                    na = 1 - na
                case (3)
                    ix2 = iw + idot
                    if (na == 0) then
                        call complex_forward_pass_3(idot, l1, c, ch, wa(iw), wa(ix2))
                    else
                        call complex_forward_pass_3(idot, l1, ch, c, wa(iw), wa(ix2))
                    end if
                    na = 1 - na
                case (4)
                    ix3 = ix2 + idot
                    if (na == 0) then
                        call complex_forward_pass_4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
                    else
                        call complex_forward_pass_4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
                    end if
                    na = 1 - na
                case (5)
                    ix2 = iw + idot
                    ix3 = ix2 + idot
                    ix4 = ix3 + idot
                    if (na == 0) then
                        call complex_forward_pass_5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    else
                        call complex_forward_pass_5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    end if
                    na = 1 - na
                case default
                    if (na == 0) then
                        call complex_forward_pass_n(nac, idot, iip, l1, idl1, c, c, c, ch, ch, wa(iw))
                    else
                        call complex_forward_pass_n(nac, idot, iip, l1, idl1, ch, ch, ch, c, c, wa(iw))
                    end if
                    if (nac /= 0) na = 1 - na
            end select
            l1 = l2
            iw = iw + (iip - 1)*idot
        end do

        if (na /= 0) then
            n2 = n + n
            c(:n2) = ch(:n2)
        end if

    end subroutine complex_forward_lower_utility_routine


    pure subroutine complex_forward_pass_2(ido, l1, cc, ch, wa1)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 2, l1)
        real(wp), intent(out) :: ch(ido, l1, 2)
        real(wp), intent(in) :: wa1(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp) :: tr2, ti2


        if(ido <= 2) then
            ch(1, :, 1) = cc(1, 1, :) + cc(1, 2, :)
            ch(1, :, 2) = cc(1, 1, :) - cc(1, 2, :)
            ch(2, :, 1) = cc(2, 1, :) + cc(2, 2, :)
            ch(2, :, 2) = cc(2, 1, :) - cc(2, 2, :)
        else
            do k = 1, l1
                do i = 2, ido, 2
                    ch(i-1, k, 1) = cc(i-1, 1, k) + cc(i-1, 2, k)
                    tr2 = cc(i-1, 1, k) - cc(i-1, 2, k)
                    ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
                    ti2 = cc(i, 1, k) - cc(i, 2, k)
                    ch(i, k, 2) = wa1(i-1)*ti2 - wa1(i)*tr2
                    ch(i-1, k, 2) = wa1(i-1)*tr2 + wa1(i)*ti2
                end do
            end do
        end if

    end subroutine complex_forward_pass_2


    pure subroutine complex_forward_pass_3(ido, l1, cc, ch, wa1, wa2)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 3, l1)
        real(wp), intent(out) :: ch(ido, l1, 3)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp), parameter :: taur = -HALF
        real(wp), parameter :: taui = -SQRT_3/2 !  - 0.866025403784439
        real(wp) :: tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3


        select case(ido)
            case (2)
                do k = 1, l1
                    tr2 = cc(1, 2, k) + cc(1, 3, k)
                    cr2 = cc(1, 1, k) + taur*tr2
                    ch(1, k, 1) = cc(1, 1, k) + tr2
                    ti2 = cc(2, 2, k) + cc(2, 3, k)
                    ci2 = cc(2, 1, k) + taur*ti2
                    ch(2, k, 1) = cc(2, 1, k) + ti2
                    cr3 = taui*(cc(1, 2, k)-cc(1, 3, k))
                    ci3 = taui*(cc(2, 2, k)-cc(2, 3, k))
                    ch(1, k, 2) = cr2 - ci3
                    ch(1, k, 3) = cr2 + ci3
                    ch(2, k, 2) = ci2 + cr3
                    ch(2, k, 3) = ci2 - cr3
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        tr2 = cc(i-1, 2, k) + cc(i-1, 3, k)
                        cr2 = cc(i-1, 1, k) + taur*tr2
                        ch(i-1, k, 1) = cc(i-1, 1, k) + tr2
                        ti2 = cc(i, 2, k) + cc(i, 3, k)
                        ci2 = cc(i, 1, k) + taur*ti2
                        ch(i, k, 1) = cc(i, 1, k) + ti2
                        cr3 = taui*(cc(i-1, 2, k)-cc(i-1, 3, k))
                        ci3 = taui*(cc(i, 2, k)-cc(i, 3, k))
                        dr2 = cr2 - ci3
                        dr3 = cr2 + ci3
                        di2 = ci2 + cr3
                        di3 = ci2 - cr3
                        ch(i, k, 2) = wa1(i-1)*di2 - wa1(i)*dr2
                        ch(i-1, k, 2) = wa1(i-1)*dr2 + wa1(i)*di2
                        ch(i, k, 3) = wa2(i-1)*di3 - wa2(i)*dr3
                        ch(i-1, k, 3) = wa2(i-1)*dr3 + wa2(i)*di3
                    end do
                end do
        end select

    end subroutine complex_forward_pass_3

    pure subroutine complex_forward_pass_4(ido, l1, cc, ch, wa1, wa2, wa3)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 4, l1)
        real(wp), intent(out) :: ch(ido, l1, 4)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp) :: ti1, ti2, tr4, ti3, tr1, tr2, ti4
        real(wp) :: tr3, cr3, ci3, cr2, cr4, ci2, ci4


        select case(ido)
            case (2)
                do k = 1, l1
                    ti1 = cc(2, 1, k) - cc(2, 3, k)
                    ti2 = cc(2, 1, k) + cc(2, 3, k)
                    tr4 = cc(2, 2, k) - cc(2, 4, k)
                    ti3 = cc(2, 2, k) + cc(2, 4, k)
                    tr1 = cc(1, 1, k) - cc(1, 3, k)
                    tr2 = cc(1, 1, k) + cc(1, 3, k)
                    ti4 = cc(1, 4, k) - cc(1, 2, k)
                    tr3 = cc(1, 2, k) + cc(1, 4, k)
                    ch(1, k, 1) = tr2 + tr3
                    ch(1, k, 3) = tr2 - tr3
                    ch(2, k, 1) = ti2 + ti3
                    ch(2, k, 3) = ti2 - ti3
                    ch(1, k, 2) = tr1 + tr4
                    ch(1, k, 4) = tr1 - tr4
                    ch(2, k, 2) = ti1 + ti4
                    ch(2, k, 4) = ti1 - ti4
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        ti1 = cc(i, 1, k) - cc(i, 3, k)
                        ti2 = cc(i, 1, k) + cc(i, 3, k)
                        ti3 = cc(i, 2, k) + cc(i, 4, k)
                        tr4 = cc(i, 2, k) - cc(i, 4, k)
                        tr1 = cc(i-1, 1, k) - cc(i-1, 3, k)
                        tr2 = cc(i-1, 1, k) + cc(i-1, 3, k)
                        ti4 = cc(i-1, 4, k) - cc(i-1, 2, k)
                        tr3 = cc(i-1, 2, k) + cc(i-1, 4, k)
                        ch(i-1, k, 1) = tr2 + tr3
                        cr3 = tr2 - tr3
                        ch(i, k, 1) = ti2 + ti3
                        ci3 = ti2 - ti3
                        cr2 = tr1 + tr4
                        cr4 = tr1 - tr4
                        ci2 = ti1 + ti4
                        ci4 = ti1 - ti4
                        ch(i-1, k, 2) = wa1(i-1)*cr2 + wa1(i)*ci2
                        ch(i, k, 2) = wa1(i-1)*ci2 - wa1(i)*cr2
                        ch(i-1, k, 3) = wa2(i-1)*cr3 + wa2(i)*ci3
                        ch(i, k, 3) = wa2(i-1)*ci3 - wa2(i)*cr3
                        ch(i-1, k, 4) = wa3(i-1)*cr4 + wa3(i)*ci4
                        ch(i, k, 4) = wa3(i-1)*ci4 - wa3(i)*cr4
                    end do
                end do
        end select

    end subroutine complex_forward_pass_4

    pure subroutine complex_forward_pass_5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 5, l1)
        real(wp), intent(out) :: ch(ido, l1, 5)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)
        real(wp), intent(in) :: wa4(*)

        ! Local variables

        integer(ip) :: k, i
        real(wp) :: ti5, ti2, ti4, ti3, tr5, tr2, tr4
        real(wp) :: tr3, cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4, dr3, dr4, di3
        real(wp) :: di4, dr5, dr2, di5, di2
        real(wp), parameter :: tr11 = (SQRT_5 - ONE)/4 ! 0.309016994374947
        real(wp), parameter :: ti11 = -sqrt((FIVE + SQRT_5)/2)/2 ! -.951056516295154
        real(wp), parameter :: tr12 = (-ONE - SQRT_5)/4 ! -.809016994374947
        real(wp), parameter :: ti12 = -sqrt(FIVE/(TWO * (FIVE + SQRT_5))) ! -0.587785252292473


        select case(ido)
            case (2)
                do k = 1, l1
                    ti5 = cc(2, 2, k) - cc(2, 5, k)
                    ti2 = cc(2, 2, k) + cc(2, 5, k)
                    ti4 = cc(2, 3, k) - cc(2, 4, k)
                    ti3 = cc(2, 3, k) + cc(2, 4, k)
                    tr5 = cc(1, 2, k) - cc(1, 5, k)
                    tr2 = cc(1, 2, k) + cc(1, 5, k)
                    tr4 = cc(1, 3, k) - cc(1, 4, k)
                    tr3 = cc(1, 3, k) + cc(1, 4, k)
                    ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
                    ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
                    cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
                    ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
                    cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
                    ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
                    cr5 = ti11*tr5 + ti12*tr4
                    ci5 = ti11*ti5 + ti12*ti4
                    cr4 = ti12*tr5 - ti11*tr4
                    ci4 = ti12*ti5 - ti11*ti4
                    ch(1, k, 2) = cr2 - ci5
                    ch(1, k, 5) = cr2 + ci5
                    ch(2, k, 2) = ci2 + cr5
                    ch(2, k, 3) = ci3 + cr4
                    ch(1, k, 3) = cr3 - ci4
                    ch(1, k, 4) = cr3 + ci4
                    ch(2, k, 4) = ci3 - cr4
                    ch(2, k, 5) = ci2 - cr5
                end do
            case default
                do k = 1, l1
                    do i = 2, ido, 2
                        ti5 = cc(i, 2, k) - cc(i, 5, k)
                        ti2 = cc(i, 2, k) + cc(i, 5, k)
                        ti4 = cc(i, 3, k) - cc(i, 4, k)
                        ti3 = cc(i, 3, k) + cc(i, 4, k)
                        tr5 = cc(i-1, 2, k) - cc(i-1, 5, k)
                        tr2 = cc(i-1, 2, k) + cc(i-1, 5, k)
                        tr4 = cc(i-1, 3, k) - cc(i-1, 4, k)
                        tr3 = cc(i-1, 3, k) + cc(i-1, 4, k)
                        ch(i-1, k, 1) = cc(i-1, 1, k) + tr2 + tr3
                        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
                        cr2 = cc(i-1, 1, k) + tr11*tr2 + tr12*tr3
                        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
                        cr3 = cc(i-1, 1, k) + tr12*tr2 + tr11*tr3
                        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
                        cr5 = ti11*tr5 + ti12*tr4
                        ci5 = ti11*ti5 + ti12*ti4
                        cr4 = ti12*tr5 - ti11*tr4
                        ci4 = ti12*ti5 - ti11*ti4
                        dr3 = cr3 - ci4
                        dr4 = cr3 + ci4
                        di3 = ci3 + cr4
                        di4 = ci3 - cr4
                        dr5 = cr2 + ci5
                        dr2 = cr2 - ci5
                        di5 = ci2 - cr5
                        di2 = ci2 + cr5
                        ch(i-1, k, 2) = wa1(i-1)*dr2 + wa1(i)*di2
                        ch(i, k, 2) = wa1(i-1)*di2 - wa1(i)*dr2
                        ch(i-1, k, 3) = wa2(i-1)*dr3 + wa2(i)*di3
                        ch(i, k, 3) = wa2(i-1)*di3 - wa2(i)*dr3
                        ch(i-1, k, 4) = wa3(i-1)*dr4 + wa3(i)*di4
                        ch(i, k, 4) = wa3(i-1)*di4 - wa3(i)*dr4
                        ch(i-1, k, 5) = wa4(i-1)*dr5 + wa4(i)*di5
                        ch(i, k, 5) = wa4(i-1)*di5 - wa4(i)*dr5
                    end do
                end do
        end select

    end subroutine complex_forward_pass_5

    subroutine complex_forward_pass_n(nac, ido, iip, l1, idl1, cc, c1, c2, ch, ch2, wa)

        ! Dummy arguments

        integer(ip), intent(out)    :: nac
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: iip
        integer(ip), intent(in)     :: l1
        integer(ip), intent(in)     :: idl1
        real(wp),    intent(inout)  :: cc(ido, iip, l1)
        real(wp),    intent(inout)  :: c1(ido, l1, iip)
        real(wp),    intent(inout)  :: c2(idl1, iip)
        real(wp),    intent(inout)  :: ch(ido, l1, iip)
        real(wp),    intent(inout)  :: ch2(idl1, iip)
        real(wp),    intent(in)     :: wa(*)

        ! Local variables

        integer(ip) :: idot, nt, iipp2, iipph, idp, j, jc
        integer(ip) :: k, i, idl, inc, l, lc
        integer(ip) :: idlj, idij, idj
        real(wp)    :: war, wai


        idot = ido/2
        nt = iip*idl1
        iipp2 = iip + 2
        iipph = (iip + 1)/2
        idp = iip*ido

        if (l1 <= ido) then
            do j = 2, iipph
                jc = iipp2 - j
                ch(:, :, j) = cc(:, j, :) + cc(:, jc, :)
                ch(:, :, jc) = cc(:, j, :) - cc(:, jc, :)
            end do
            ch(:, :, 1) = cc(:, 1, :)
        else
            do j = 2, iipph
                jc = iipp2 - j
                ch(:, :, j) = cc(:, j, :) + cc(:, jc, :)
                ch(:, :, jc) = cc(:, j, :) - cc(:, jc, :)
            end do
            ch(:, :, 1) = cc(:, 1, :)
        end if

        idl = 2 - ido
        inc = 0

        do l = 2, iipph
            lc = iipp2 - l
            idl = idl + ido
            c2(:, l) = ch2(:, 1) + wa(idl-1)*ch2(:, 2)
            c2(:, lc) = -wa(idl)*ch2(:, iip)
            idlj = idl
            inc = inc + ido
            do j = 3, iipph
                jc = iipp2 - j
                idlj = idlj + inc
                if (idlj > idp) then
                    idlj = idlj - idp
                end if
                war = wa(idlj-1)
                wai = wa(idlj)
                c2(:, l) = c2(:, l) + war*ch2(:, j)
                c2(:, lc) = c2(:, lc) - wai*ch2(:, jc)
            end do
        end do

        do j = 2, iipph
            ch2(:, 1) = ch2(:, 1) + ch2(:, j)
        end do

        do j = 2, iipph
            jc = iipp2 - j
            ch2(:idl1-1:2, j) = c2(:idl1-1:2, j) - c2(2:idl1:2, jc)
            ch2(:idl1-1:2, jc) = c2(:idl1-1:2, j) + c2(2:idl1:2, jc)
            ch2(2:idl1:2, j) = c2(2:idl1:2, j) + c2(:idl1-1:2, jc)
            ch2(2:idl1:2, jc) = c2(2:idl1:2, j) - c2(:idl1-1:2, jc)
        end do

        nac = 1

        if (ido /= 2) then

            nac = 0
            c2(:, 1) = ch2(:, 1)
            c1(1, :, 2:iip) = ch(1, :, 2:iip)
            c1(2, :, 2:iip) = ch(2, :, 2:iip)

            if (idot <= l1) then
                idij = 0
                do j = 2, iip
                    idij = idij + 2
                    do i = 4, ido, 2
                        idij = idij + 2
                        c1(i-1, :, j) = wa(idij-1)*ch(i-1, :, j) + wa(idij)*ch(i, :, j)
                        c1(i, :, j) = wa(idij-1)*ch(i, :, j) - wa(idij)*ch(i-1, :, j)
                    end do
                end do
            else
                idj = 2 - ido
                do j = 2, iip
                    idj = idj + ido
                    do k = 1, l1
                        idij = idj
                        c1(3:ido-1:2, k, j) = wa(idij+1:ido-3+idij:2)*ch(3:ido-1:2, k, j &
                            ) + wa(idij+2:ido-2+idij:2)*ch(4:ido:2, k, j)
                        c1(4:ido:2, k, j) = wa(idij+1:ido-3+idij:2)*ch(4:ido:2, k, j) - &
                            wa(idij+2:ido-2+idij:2)*ch(3:ido-1:2, k, j)
                    end do
                end do
            end if
        end if

    end subroutine complex_forward_pass_n

    subroutine rffti(n, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: wsave(*)


        if (n /= 1) then
            call real_initialization_lower_utility_routine(n, wsave(n+1), wsave(2*n+1))
        end if

    end subroutine rffti

    subroutine real_initialization_lower_utility_routine(n, wa, ifac)

        ! Dummy arguments

        integer(ip),  intent(in)     :: n
        real(wp), intent(inout)  :: ifac(*)
        real(wp), intent(out)    :: wa(*)

        ! Local variables

        integer, parameter :: NTRYH(*) = [4, 2, 3, 5]
        integer(ip) :: nl, nf, j, ntry, nq, nr, i, is, nfm1, l1, k1, iip
        integer(ip) :: ld, l2, ido, iipm, ii
        real(wp) :: argh, argld, fi, arg


        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorization_loop: do

            j = j + 1

            if (j <= 4) then
                ntry = NTRYH(j)
            else
                ntry = ntry + 2
            end if

            inner_loop: do

                nq = nl/ntry
                nr = nl - ntry*nq

                if (nr /= 0) cycle factorization_loop

                nf = nf + 1
                ifac(nf+2) = ntry
                nl = nq

                if (ntry == 2) then
                    if (nf /= 1) then
                        ifac(nf+2:4:(-1)) = ifac(nf+1:3:(-1))
                        ifac(3) = 2
                    end if
                end if

                if (nl /= 1) cycle inner_loop

                exit inner_loop
            end do inner_loop
            exit factorization_loop
        end do factorization_loop

        ifac(1) = n
        ifac(2) = nf
        argh = TWO_PI/n
        is = 0
        nfm1 = nf - 1
        l1 = 1

        if (nfm1 /= 0) then
            do k1 = 1, nfm1
                iip = int(ifac(k1+2), kind=ip)
                ld = 0
                l2 = l1*iip
                ido = n/l2
                iipm = iip - 1
                do j = 1, iipm
                    ld = ld + l1
                    i = is
                    argld = real(ld, kind=wp)*argh
                    fi = ZERO
                    do ii = 3, ido, 2
                        i = i + 2
                        fi = fi + ONE
                        arg = fi*argld
                        wa(i-1) = cos(arg)
                        wa(i) = sin(arg)
                    end do
                    is = is + ido
                end do
                l1 = l2
            end do
        end if

    end subroutine real_initialization_lower_utility_routine

    subroutine rfftb(n, r, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: r(*)
        real(wp) :: wsave(*)


        if (n /= 1) then
            call real_backward_lower_utility_routine(n, r, wsave, wsave(n+1), wsave(2*n+1))
        end if

    end subroutine rfftb

    subroutine real_backward_lower_utility_routine(n, c, ch, wa, ifac)

        ! Dummy arguments

        integer(ip),  intent(in) :: n
        real(wp), intent(in) :: ifac(*)
        real(wp)              :: c(*)
        real(wp)              :: ch(*)
        real(wp)              :: wa(*)

        ! Local variables

        integer(ip) :: nf, na, l1, iw, k1, iip
        integer(ip) :: l2, ido, idl1, ix2, ix3, ix4


        nf = int(ifac(2), kind=ip)
        na = 0
        l1 = 1
        iw = 1

        do k1 = 1, nf
            iip = int(ifac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idl1 = ido*l1

            select case (iip)
                case (4)
                    ix2 = iw + ido
                    ix3 = ix2 + ido
                    if (na == 0) then
                        call real_backward_pass_4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
                    else
                        call real_backward_pass_4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
                    end if
                    na = 1 - na
                case (2)
                    if (na == 0) then
                        call real_backward_pass_2(ido, l1, c, ch, wa(iw))
                    else
                        call real_backward_pass_2(ido, l1, ch, c, wa(iw))
                    end if
                    na = 1 - na
                case (3)
                    ix2 = iw + ido
                    if (na == 0) then
                        call real_backward_pass_3(ido, l1, c, ch, wa(iw), wa(ix2))
                    else
                        call real_backward_pass_3(ido, l1, ch, c, wa(iw), wa(ix2))
                    end if
                    na = 1 - na
                case (5)
                    ix2 = iw + ido
                    ix3 = ix2 + ido
                    ix4 = ix3 + ido
                    if (na == 0) then
                        call real_backward_pass_5(ido, l1, c, ch, wa(iw), wa(ix2), wa( &
                            ix3), wa(ix4))
                    else
                        call real_backward_pass_5(ido, l1, ch, c, wa(iw), wa(ix2), wa( &
                            ix3), wa(ix4))
                    end if
                    na = 1 - na
                case default
                    if (na == 0) then
                        call real_backward_pass_n(ido, iip, l1, idl1, c, c, c, ch, ch, wa(iw))
                    else
                        call real_backward_pass_n(ido, iip, l1, idl1, ch, ch, ch, c, c, wa(iw))
                    end if

                    if (ido == 1) na = 1 - na

            end select
            l1 = l2
            iw = iw + (iip - 1)*ido
        end do

        if (na /= 0) c(:n) = ch(:n)

    end subroutine real_backward_lower_utility_routine

    pure subroutine real_backward_pass_2(ido, l1, cc, ch, wa1)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 2, l1)
        real(wp), intent(out) :: ch(ido, l1, 2)
        real(wp), intent(in) :: wa1(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp) :: tr2, ti2


        ch(1, :, 1) = cc(1, 1, :) + cc(ido, 2, :)
        ch(1, :, 2) = cc(1, 1, :) - cc(ido, 2, :)

        if (ido >= 2) then
            if (ido /= 2) then
                idp2 = ido + 2
                do k = 1, l1
                    do i = 3, ido, 2
                        ic = idp2 - i
                        ch(i-1, k, 1) = cc(i-1, 1, k) + cc(ic-1, 2, k)
                        tr2 = cc(i-1, 1, k) - cc(ic-1, 2, k)
                        ch(i, k, 1) = cc(i, 1, k) - cc(ic, 2, k)
                        ti2 = cc(i, 1, k) + cc(ic, 2, k)
                        ch(i-1, k, 2) = wa1(i-2)*tr2 - wa1(i-1)*ti2
                        ch(i, k, 2) = wa1(i-2)*ti2 + wa1(i-1)*tr2
                    end do
                end do
                if (mod(ido, 2) == 1) then
                    return
                end if
            end if
            ch(ido, :, 1) = cc(ido, 1, :) + cc(ido, 1, :)
            ch(ido, :, 2) = -(cc(1, 2, :)+cc(1, 2, :))
        end if

    end subroutine real_backward_pass_2

    pure subroutine real_backward_pass_3(ido, l1, cc, ch, wa1, wa2)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 3, l1)
        real(wp), intent(out) :: ch(ido, l1, 3)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp), parameter :: taur = -HALF
        real(wp), parameter :: taui = SQRT_3/TWO ! 0.866025403784439
        real(wp)            :: tr2, cr2, ci3, ti2, ci2, cr3, dr2, dr3, di2, di3


        do k = 1, l1
            tr2 = cc(ido, 2, k) + cc(ido, 2, k)
            cr2 = cc(1, 1, k) + taur*tr2
            ch(1, k, 1) = cc(1, 1, k) + tr2
            ci3 = taui*(cc(1, 3, k)+cc(1, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
        end do

        if(ido /= 1) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    tr2 = cc(i-1, 3, k) + cc(ic-1, 2, k)
                    cr2 = cc(i-1, 1, k) + taur*tr2
                    ch(i-1, k, 1) = cc(i-1, 1, k) + tr2
                    ti2 = cc(i, 3, k) - cc(ic, 2, k)
                    ci2 = cc(i, 1, k) + taur*ti2
                    ch(i, k, 1) = cc(i, 1, k) + ti2
                    cr3 = taui*(cc(i-1, 3, k)-cc(ic-1, 2, k))
                    ci3 = taui*(cc(i, 3, k)+cc(ic, 2, k))
                    dr2 = cr2 - ci3
                    dr3 = cr2 + ci3
                    di2 = ci2 + cr3
                    di3 = ci2 - cr3
                    ch(i-1, k, 2) = wa1(i-2)*dr2 - wa1(i-1)*di2
                    ch(i, k, 2) = wa1(i-2)*di2 + wa1(i-1)*dr2
                    ch(i-1, k, 3) = wa2(i-2)*dr3 - wa2(i-1)*di3
                    ch(i, k, 3) = wa2(i-2)*di3 + wa2(i-1)*dr3
                end do
            end do
        end if

    end subroutine real_backward_pass_3

    pure subroutine real_backward_pass_4(ido, l1, cc, ch, wa1, wa2, wa3)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 4, l1)
        real(wp), intent(out) :: ch(ido, l1, 4)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)

        ! Local variables

        integer(ip)         :: k, idp2, i, ic
        real(wp), parameter :: SQRT_2 = sqrt(TWO) ! 1.414213562373095
        real(wp)            :: tr1, tr2, tr3, tr4, ti1, ti2, ti3, ti4, cr3, ci3
        real(wp)            :: cr2, cr4, ci2, ci4


        do k = 1, l1
            tr1 = cc(1, 1, k) - cc(ido, 4, k)
            tr2 = cc(1, 1, k) + cc(ido, 4, k)
            tr3 = cc(ido, 2, k) + cc(ido, 2, k)
            tr4 = cc(1, 3, k) + cc(1, 3, k)
            ch(1, k, 1) = tr2 + tr3
            ch(1, k, 2) = tr1 - tr4
            ch(1, k, 3) = tr2 - tr3
            ch(1, k, 4) = tr1 + tr4
        end do

        if (2 <= ido) then
            if (ido /= 2) then
                idp2 = ido + 2
                do k = 1, l1
                    do i = 3, ido, 2
                        ic = idp2 - i
                        ti1 = cc(i, 1, k) + cc(ic, 4, k)
                        ti2 = cc(i, 1, k) - cc(ic, 4, k)
                        ti3 = cc(i, 3, k) - cc(ic, 2, k)
                        tr4 = cc(i, 3, k) + cc(ic, 2, k)
                        tr1 = cc(i-1, 1, k) - cc(ic-1, 4, k)
                        tr2 = cc(i-1, 1, k) + cc(ic-1, 4, k)
                        ti4 = cc(i-1, 3, k) - cc(ic-1, 2, k)
                        tr3 = cc(i-1, 3, k) + cc(ic-1, 2, k)
                        ch(i-1, k, 1) = tr2 + tr3
                        cr3 = tr2 - tr3
                        ch(i, k, 1) = ti2 + ti3
                        ci3 = ti2 - ti3
                        cr2 = tr1 - tr4
                        cr4 = tr1 + tr4
                        ci2 = ti1 + ti4
                        ci4 = ti1 - ti4
                        ch(i-1, k, 2) = wa1(i-2)*cr2 - wa1(i-1)*ci2
                        ch(i, k, 2) = wa1(i-2)*ci2 + wa1(i-1)*cr2
                        ch(i-1, k, 3) = wa2(i-2)*cr3 - wa2(i-1)*ci3
                        ch(i, k, 3) = wa2(i-2)*ci3 + wa2(i-1)*cr3
                        ch(i-1, k, 4) = wa3(i-2)*cr4 - wa3(i-1)*ci4
                        ch(i, k, 4) = wa3(i-2)*ci4 + wa3(i-1)*cr4
                    end do
                end do
                if (mod(ido, 2) == 1) then
                    return
                end if
            end if
            do k = 1, l1
                ti1 = cc(1, 2, k) + cc(1, 4, k)
                ti2 = cc(1, 4, k) - cc(1, 2, k)
                tr1 = cc(ido, 1, k) - cc(ido, 3, k)
                tr2 = cc(ido, 1, k) + cc(ido, 3, k)
                ch(ido, k, 1) = tr2 + tr2
                ch(ido, k, 2) = SQRT_2*(tr1 - ti1)
                ch(ido, k, 3) = ti2 + ti2
                ch(ido, k, 4) = -SQRT_2*(tr1 + ti1)
            end do
        end if

    end subroutine real_backward_pass_4

    pure subroutine real_backward_pass_5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, 5, l1)
        real(wp), intent(out) :: ch(ido, l1, 5)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)
        real(wp), intent(in) :: wa3(*)
        real(wp), intent(in) :: wa4(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp) :: ti5, ti4, tr2, tr3, cr2, cr3, ci5
        real(wp) :: ci4, ti2, ti3, tr5, tr4, ci2, ci3, cr5, cr4, dr3, dr4, di3
        real(wp) :: di4, dr5, dr2, di5, di2
        real(wp), parameter :: tr11 = (SQRT_5 - ONE)/4 ! 0.309016994374947
        real(wp), parameter :: ti11 = sqrt((FIVE + SQRT_5)/2)/2 ! 0.951056516295154
        real(wp), parameter :: tr12 = (-ONE - SQRT_5)/4 ! -.809016994374947
        real(wp), parameter :: ti12 = sqrt(FIVE/(TWO * (FIVE + SQRT_5))) ! 0.587785252292473



        do k = 1, l1
            ti5 = cc(1, 3, k) + cc(1, 3, k)
            ti4 = cc(1, 5, k) + cc(1, 5, k)
            tr2 = cc(ido, 2, k) + cc(ido, 2, k)
            tr3 = cc(ido, 4, k) + cc(ido, 4, k)
            ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
            cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
            cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
            ci5 = ti11*ti5 + ti12*ti4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1, k, 2) = cr2 - ci5
            ch(1, k, 3) = cr3 - ci4
            ch(1, k, 4) = cr3 + ci4
            ch(1, k, 5) = cr2 + ci5
        end do

        if(ido /= 1) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    ti5 = cc(i, 3, k) + cc(ic, 2, k)
                    ti2 = cc(i, 3, k) - cc(ic, 2, k)
                    ti4 = cc(i, 5, k) + cc(ic, 4, k)
                    ti3 = cc(i, 5, k) - cc(ic, 4, k)
                    tr5 = cc(i-1, 3, k) - cc(ic-1, 2, k)
                    tr2 = cc(i-1, 3, k) + cc(ic-1, 2, k)
                    tr4 = cc(i-1, 5, k) - cc(ic-1, 4, k)
                    tr3 = cc(i-1, 5, k) + cc(ic-1, 4, k)
                    ch(i-1, k, 1) = cc(i-1, 1, k) + tr2 + tr3
                    ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
                    cr2 = cc(i-1, 1, k) + tr11*tr2 + tr12*tr3
                    ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
                    cr3 = cc(i-1, 1, k) + tr12*tr2 + tr11*tr3
                    ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
                    cr5 = ti11*tr5 + ti12*tr4
                    ci5 = ti11*ti5 + ti12*ti4
                    cr4 = ti12*tr5 - ti11*tr4
                    ci4 = ti12*ti5 - ti11*ti4
                    dr3 = cr3 - ci4
                    dr4 = cr3 + ci4
                    di3 = ci3 + cr4
                    di4 = ci3 - cr4
                    dr5 = cr2 + ci5
                    dr2 = cr2 - ci5
                    di5 = ci2 - cr5
                    di2 = ci2 + cr5
                    ch(i-1, k, 2) = wa1(i-2)*dr2 - wa1(i-1)*di2
                    ch(i, k, 2) = wa1(i-2)*di2 + wa1(i-1)*dr2
                    ch(i-1, k, 3) = wa2(i-2)*dr3 - wa2(i-1)*di3
                    ch(i, k, 3) = wa2(i-2)*di3 + wa2(i-1)*dr3
                    ch(i-1, k, 4) = wa3(i-2)*dr4 - wa3(i-1)*di4
                    ch(i, k, 4) = wa3(i-2)*di4 + wa3(i-1)*dr4
                    ch(i-1, k, 5) = wa4(i-2)*dr5 - wa4(i-1)*di5
                    ch(i, k, 5) = wa4(i-2)*di5 + wa4(i-1)*dr5
                end do
            end do
        end if

    end subroutine real_backward_pass_5

    subroutine real_backward_pass_n(ido, iip, l1, idl1, cc, c1, c2, ch, ch2, wa)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: iip
        integer(ip), intent(in) :: l1
        integer(ip), intent(in) :: idl1
        real(wp), intent(in) :: cc(ido, iip, l1)
        real(wp), intent(inout)  :: c1(ido, l1, iip)
        real(wp), intent(inout)  :: c2(idl1, iip)
        real(wp), intent(inout)  :: ch(ido, l1, iip)
        real(wp), intent(inout)  :: ch2(idl1, iip)
        real(wp), intent(in) :: wa(*)

        ! Local variables

        integer:: idp2, nbd, iipp2, iipph, k, i, j, jc, j2, l, lc, is, idij
        real(wp) :: arg, dcp, dsp, ar1, ai1, ar1h, dc2, ds2, ar2, ai2, ar2h


        arg = TWO_PI/iip
        dcp = cos(arg)
        dsp = sin(arg)
        idp2 = ido + 2
        nbd =(ido - 1)/2
        iipp2 = iip + 2
        iipph = (iip + 1)/2

        ch(:, :, 1) = cc(:, 1, :)

        do j = 2, iipph
            jc = iipp2 - j
            j2 = j + j
            ch(1, :, j) = cc(ido, j2-2, :) + cc(ido, j2-2, :)
            ch(1, :, jc) = cc(1, j2-1, :) + cc(1, j2-1, :)
        end do

        if(ido /= 1) then
            if (l1 <= nbd) then
                do j = 2, iipph
                    jc = iipp2 - j
                    ch(2:ido-1:2, :, j) = cc(2:ido-1:2, 2*j-1, :) + cc(idp2-4: &
                        idp2-1-ido:(-2), 2*j-2, :)
                    ch(2:ido-1:2, :, jc) = cc(2:ido-1:2, 2*j-1, :) - cc(idp2-4: &
                        idp2-1-ido:(-2), 2*j-2, :)
                    ch(3:ido:2, :, j) = cc(3:ido:2, 2*j-1, :) - cc(idp2-3:idp2- &
                        ido:(-2), 2*j-2, :)
                    ch(3:ido:2, :, jc) = cc(3:ido:2, 2*j-1, :) + cc(idp2-3:idp2- &
                        ido:(-2), 2*j-2, :)
                end do
            else
                do j = 2, iipph
                    jc = iipp2 - j
                    ch(2:ido-1:2, :, j) = cc(2:ido-1:2, 2*j-1, :) + cc(idp2-4: &
                        idp2-1-ido:(-2), 2*j-2, :)
                    ch(2:ido-1:2, :, jc) = cc(2:ido-1:2, 2*j-1, :) - cc(idp2-4: &
                        idp2-1-ido:(-2), 2*j-2, :)
                    ch(3:ido:2, :, j) = cc(3:ido:2, 2*j-1, :) - cc(idp2-3:idp2- &
                        ido:(-2), 2*j-2, :)
                    ch(3:ido:2, :, jc) = cc(3:ido:2, 2*j-1, :) + cc(idp2-3:idp2- &
                        ido:(-2), 2*j-2, :)
                end do
            end if
        end if

        ar1 = ONE
        ai1 = ZERO

        do l = 2, iipph
            lc = iipp2 - l
            ar1h = dcp*ar1 - dsp*ai1
            ai1 = dcp*ai1 + dsp*ar1
            ar1 = ar1h
            c2(:, l) = ch2(:, 1) + ar1*ch2(:, 2)
            c2(:, lc) = ai1*ch2(:, iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j = 3, iipph
                jc = iipp2 - j
                ar2h = dc2*ar2 - ds2*ai2
                ai2 = dc2*ai2 + ds2*ar2
                ar2 = ar2h
                c2(:, l) = c2(:, l) + ar2*ch2(:, j)
                c2(:, lc) = c2(:, lc) + ai2*ch2(:, jc)
            end do
        end do

        do j = 2, iipph
            ch2(:, 1) = ch2(:, 1) + ch2(:, j)
        end do

        do j = 2, iipph
            jc = iipp2 - j
            ch(1, :, j) = c1(1, :, j) - c1(1, :, jc)
            ch(1, :, jc) = c1(1, :, j) + c1(1, :, jc)
        end do

        if(ido /= 1) then
            if (l1 <= nbd) then
                do j = 2, iipph
                    jc = iipp2 - j
                    ch(2:ido-1:2, :, j) = c1(2:ido-1:2, :, j) - c1(3:ido:2, :, jc)
                    ch(2:ido-1:2, :, jc) = c1(2:ido-1:2, :, j) + c1(3:ido:2, :, jc)
                    ch(3:ido:2, :, j) = c1(3:ido:2, :, j) + c1(2:ido-1:2, :, jc)
                    ch(3:ido:2, :, jc) = c1(3:ido:2, :, j) - c1(2:ido-1:2, :, jc)
                end do
            else
                do j = 2, iipph
                    jc = iipp2 - j
                    ch(2:ido-1:2, :, j) = c1(2:ido-1:2, :, j) - c1(3:ido:2, :, jc)
                    ch(2:ido-1:2, :, jc) = c1(2:ido-1:2, :, j) + c1(3:ido:2, :, jc)
                    ch(3:ido:2, :, j) = c1(3:ido:2, :, j) + c1(2:ido-1:2, :, jc)
                    ch(3:ido:2, :, jc) = c1(3:ido:2, :, j) - c1(2:ido-1:2, :, jc)
                end do
            end if
        end if

        if(ido /= 1) then
            c2(:, 1) = ch2(:, 1)
            c1(1, :, 2:iip) = ch(1, :, 2:iip)
            if (nbd <= l1) then
                is = -ido
                do j = 2, iip
                    is = is + ido
                    idij = is
                    do i = 3, ido, 2
                        idij = idij + 2
                        c1(i-1, :, j) = wa(idij-1)*ch(i-1, :, j) - wa(idij)*ch(i, :, j)
                        c1(i, :, j) = wa(idij-1)*ch(i, :, j) + wa(idij)*ch(i-1, :, j)
                    end do
                end do
            else
                is = -ido
                do j = 2, iip
                    is = is + ido
                    do k = 1, l1
                        idij = is
                        c1(2:ido-1:2, k, j) = wa(idij+1:ido-2+idij:2)*ch(2:ido-1:2, &
                            k, j) - wa(idij+2:ido-1+idij:2)*ch(3:ido:2, k, j)
                        c1(3:ido:2, k, j) = wa(idij+1:ido-2+idij:2)*ch(3:ido:2, k, j) &
                            + wa(idij+2:ido-1+idij:2)*ch(2:ido-1:2, k, j)
                    end do
                end do
            end if
        end if

    end subroutine real_backward_pass_n

    subroutine rfftf(n, r, wsave)

        ! Dummy arguments

        integer(ip) :: n
        real(wp) :: r(*)
        real(wp) :: wsave(*)


        if (n /= 1) then
            call real_forward_lower_utility_routine(n, r, wsave, wsave(n+1), wsave(2*n+1))
        end if

    end subroutine rfftf

    subroutine real_forward_lower_utility_routine(n, c, ch, wa, ifac)

        ! Dummy arguments

        integer(ip),  intent(in) :: n
        real(wp), intent(in) :: ifac(*)
        real(wp)              :: c(*)
        real(wp)              :: ch(*)
        real(wp)              :: wa(*)

        ! Local variables

        integer(ip) :: nf, na, l2, iw, k1, kh
        integer(ip) :: iip, l1, ido, idl1, ix2, ix3, ix4


        nf = int(ifac(2), kind=ip)
        na = 1
        l2 = n
        iw = n

        do k1 = 1, nf
            kh = nf - k1
            iip = int(ifac(kh+3), kind=ip)
            l1 = l2/iip
            ido = n/l2
            idl1 = ido*l1
            iw = iw - (iip - 1)*ido
            na = 1 - na

            select case (iip)
                case (2)
                    if (na == 0) then
                        call real_forward_pass_2(ido, l1, c, ch, wa(iw))
                    else
                        call real_forward_pass_2(ido, l1, ch, c, wa(iw))
                    end if
                case (3)
                    ix2 = iw + ido
                    if (na == 0) then
                        call real_forward_pass_3(ido, l1, c, ch, wa(iw), wa(ix2))
                    else
                        call real_forward_pass_3(ido, l1, ch, c, wa(iw), wa(ix2))
                    end if
                case (4)
                    ix2 = iw + ido
                    ix3 = ix2 + ido
                    if (na == 0) then
                        call real_forward_pass_4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
                    else
                        call real_forward_pass_4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
                    end if
                case (5)
                    ix2 = iw + ido
                    ix3 = ix2 + ido
                    ix4 = ix3 + ido
                    if (na == 0) then
                        call real_forward_pass_5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    else
                        call real_forward_pass_5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    end if
                case default
                    if(ido == 1) then
                        na = 1 - na
                    end if
                    if (na == 0) then
                        call real_forward_pass_n(ido, iip, l1, idl1, c, c, c, ch, ch, wa(iw))
                        na = 1
                    else
                        call real_forward_pass_n(ido, iip, l1, idl1, ch, ch, ch, c, c, wa(iw))
                        na = 0
                    end if
            end select
            l2 = l1
        end do

        if (na /= 1) c(:n) = ch(:n)

    end subroutine real_forward_lower_utility_routine

    pure subroutine real_forward_pass_2(ido, l1, cc, ch, wa1)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, l1, 2)
        real(wp), intent(out) :: ch(ido, 2, l1)
        real(wp), intent(in) :: wa1(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp) :: tr2, ti2


        ch(1, 1, :) = cc(1, :, 1) + cc(1, :, 2)
        ch(ido, 2, :) = cc(1, :, 1) - cc(1, :, 2)

        if (ido >= 2) then
            if (ido /= 2) then
                idp2 = ido + 2
                do k = 1, l1
                    do i = 3, ido, 2
                        ic = idp2 - i
                        tr2 = wa1(i-2)*cc(i-1, k, 2) + wa1(i-1)*cc(i, k, 2)
                        ti2 = wa1(i-2)*cc(i, k, 2) - wa1(i-1)*cc(i-1, k, 2)
                        ch(i, 1, k) = cc(i, k, 1) + ti2
                        ch(ic, 2, k) = ti2 - cc(i, k, 1)
                        ch(i-1, 1, k) = cc(i-1, k, 1) + tr2
                        ch(ic-1, 2, k) = cc(i-1, k, 1) - tr2
                    end do
                end do
                if (mod(ido, 2) == 1) then
                    return
                end if
            end if
            ch(1, 2, :) = -cc(ido, :, 2)
            ch(ido, 1, :) = cc(ido, :, 1)
        end if

    end subroutine real_forward_pass_2

    pure subroutine real_forward_pass_3(ido, l1, cc, ch, wa1, wa2)

        ! Dummy arguments

        integer(ip), intent(in) :: ido
        integer(ip), intent(in) :: l1
        real(wp), intent(in) :: cc(ido, l1, 3)
        real(wp), intent(out) :: ch(ido, 3, l1)
        real(wp), intent(in) :: wa1(*)
        real(wp), intent(in) :: wa2(*)

        ! Local variables

        integer(ip)         :: k, idp2, i, ic
        real(wp), parameter :: taur = -HALF
        real(wp), parameter :: taui = SQRT_3/TWO ! 0.866025403784439
        real(wp)            :: cr2, dr2, di2, dr3, di3, ci2, tr2, ti2, tr3, ti3


        do k = 1, l1
            cr2 = cc(1, k, 2) + cc(1, k, 3)
            ch(1, 1, k) = cc(1, k, 1) + cr2
            ch(1, 3, k) = taui*(cc(1, k, 3)-cc(1, k, 2))
            ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
        end do

        if(ido /= 1) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    dr2 = wa1(i-2)*cc(i-1, k, 2) + wa1(i-1)*cc(i, k, 2)
                    di2 = wa1(i-2)*cc(i, k, 2) - wa1(i-1)*cc(i-1, k, 2)
                    dr3 = wa2(i-2)*cc(i-1, k, 3) + wa2(i-1)*cc(i, k, 3)
                    di3 = wa2(i-2)*cc(i, k, 3) - wa2(i-1)*cc(i-1, k, 3)
                    cr2 = dr2 + dr3
                    ci2 = di2 + di3
                    ch(i-1, 1, k) = cc(i-1, k, 1) + cr2
                    ch(i, 1, k) = cc(i, k, 1) + ci2
                    tr2 = cc(i-1, k, 1) + taur*cr2
                    ti2 = cc(i, k, 1) + taur*ci2
                    tr3 = taui*(di2 - di3)
                    ti3 = taui*(dr3 - dr2)
                    ch(i-1, 3, k) = tr2 + tr3
                    ch(ic-1, 2, k) = tr2 - tr3
                    ch(i, 3, k) = ti2 + ti3
                    ch(ic, 2, k) = ti3 - ti2
                end do
            end do
        end if

    end subroutine real_forward_pass_3

    subroutine real_forward_pass_4(ido, l1, cc, ch, wa1, wa2, wa3)

        ! Dummy arguments

        integer(ip),  intent(in)  :: ido
        integer(ip),  intent(in)  :: l1
        real(wp), intent(in)  :: cc(ido, l1, 4)
        real(wp), intent(out) :: ch(ido, 4, l1)
        real(wp), intent(in)  :: wa1(*)
        real(wp), intent(in)  :: wa2(*)
        real(wp), intent(in)  :: wa3(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp) :: tr1, tr2, cr2, ci2, cr3, ci3
        real(wp) :: cr4, ci4, tr4, ti1
        real(wp) :: ti4, ti2, ti3, tr3
        real(wp), parameter :: ONE_OVER_SQRT_2 = ONE/SQRT_2


        do k = 1, l1
            tr1 = cc(1, k, 2) + cc(1, k, 4)
            tr2 = cc(1, k, 1) + cc(1, k, 3)
            ch(1, 1, k) = tr1 + tr2
            ch(ido, 4, k) = tr2 - tr1
            ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
            ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
        end do

        if (ido >= 2) then
            if (ido /= 2) then
                idp2 = ido + 2
                do k = 1, l1
                    do i = 3, ido, 2
                        ic = idp2 - i
                        cr2 = wa1(i-2)*cc(i-1, k, 2) + wa1(i-1)*cc(i, k, 2)
                        ci2 = wa1(i-2)*cc(i, k, 2) - wa1(i-1)*cc(i-1, k, 2)
                        cr3 = wa2(i-2)*cc(i-1, k, 3) + wa2(i-1)*cc(i, k, 3)
                        ci3 = wa2(i-2)*cc(i, k, 3) - wa2(i-1)*cc(i-1, k, 3)
                        cr4 = wa3(i-2)*cc(i-1, k, 4) + wa3(i-1)*cc(i, k, 4)
                        ci4 = wa3(i-2)*cc(i, k, 4) - wa3(i-1)*cc(i-1, k, 4)
                        tr1 = cr2 + cr4
                        tr4 = cr4 - cr2
                        ti1 = ci2 + ci4
                        ti4 = ci2 - ci4
                        ti2 = cc(i, k, 1) + ci3
                        ti3 = cc(i, k, 1) - ci3
                        tr2 = cc(i-1, k, 1) + cr3
                        tr3 = cc(i-1, k, 1) - cr3
                        ch(i-1, 1, k) = tr1 + tr2
                        ch(ic-1, 4, k) = tr2 - tr1
                        ch(i, 1, k) = ti1 + ti2
                        ch(ic, 4, k) = ti1 - ti2
                        ch(i-1, 3, k) = ti4 + tr3
                        ch(ic-1, 2, k) = tr3 - ti4
                        ch(i, 3, k) = tr4 + ti3
                        ch(ic, 2, k) = tr4 - ti3
                    end do
                end do
                if (mod(ido, 2) == 1) then
                    return
                end if
            end if
            do k = 1, l1
                ti1 = -ONE_OVER_SQRT_2*(cc(ido, k, 2)+cc(ido, k, 4))
                tr1 = ONE_OVER_SQRT_2*(cc(ido, k, 2)-cc(ido, k, 4))
                ch(ido, 1, k) = tr1 + cc(ido, k, 1)
                ch(ido, 3, k) = cc(ido, k, 1) - tr1
                ch(1, 2, k) = ti1 - cc(ido, k, 3)
                ch(1, 4, k) = ti1 + cc(ido, k, 3)
            end do
        end if

    end subroutine real_forward_pass_4

    subroutine real_forward_pass_5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)

        ! Dummy arguments

        integer(ip),  intent(in)  :: ido
        integer(ip),  intent(in)  :: l1
        real(wp), intent(in)  :: cc(ido, l1, 5)
        real(wp), intent(out) :: ch(ido, 5, l1)
        real(wp), intent(in)  :: wa1(*)
        real(wp), intent(in)  :: wa2(*)
        real(wp), intent(in)  :: wa3(*)
        real(wp), intent(in)  :: wa4(*)

        ! Local variables

        integer(ip) :: k, idp2, i, ic
        real(wp) :: cr2, ci5, cr3, ci4, dr2, di2, dr3
        real(wp) :: di3, dr4, di4, dr5, di5, cr5, ci2, cr4, ci3, tr2, ti2, tr3
        real(wp) :: ti3, tr5, ti5, tr4, ti4
        real(wp), parameter :: tr11 = (SQRT_5 - ONE)/4
        real(wp), parameter :: ti11 = sqrt((FIVE + SQRT_5)/2)/2
        real(wp), parameter :: tr12 = (-ONE - SQRT_5)/4
        real(wp), parameter :: ti12 = sqrt(FIVE/(TWO * (FIVE + SQRT_5)))


        do k = 1, l1
            cr2 = cc(1, k, 5) + cc(1, k, 2)
            ci5 = cc(1, k, 5) - cc(1, k, 2)
            cr3 = cc(1, k, 4) + cc(1, k, 3)
            ci4 = cc(1, k, 4) - cc(1, k, 3)
            ch(1, 1, k) = cc(1, k, 1) + cr2 + cr3
            ch(ido, 2, k) = cc(1, k, 1) + tr11*cr2 + tr12*cr3
            ch(1, 3, k) = ti11*ci5 + ti12*ci4
            ch(ido, 4, k) = cc(1, k, 1) + tr12*cr2 + tr11*cr3
            ch(1, 5, k) = ti12*ci5 - ti11*ci4
        end do

        if(ido /= 1) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    dr2 = wa1(i-2)*cc(i-1, k, 2) + wa1(i-1)*cc(i, k, 2)
                    di2 = wa1(i-2)*cc(i, k, 2) - wa1(i-1)*cc(i-1, k, 2)
                    dr3 = wa2(i-2)*cc(i-1, k, 3) + wa2(i-1)*cc(i, k, 3)
                    di3 = wa2(i-2)*cc(i, k, 3) - wa2(i-1)*cc(i-1, k, 3)
                    dr4 = wa3(i-2)*cc(i-1, k, 4) + wa3(i-1)*cc(i, k, 4)
                    di4 = wa3(i-2)*cc(i, k, 4) - wa3(i-1)*cc(i-1, k, 4)
                    dr5 = wa4(i-2)*cc(i-1, k, 5) + wa4(i-1)*cc(i, k, 5)
                    di5 = wa4(i-2)*cc(i, k, 5) - wa4(i-1)*cc(i-1, k, 5)
                    cr2 = dr2 + dr5
                    ci5 = dr5 - dr2
                    cr5 = di2 - di5
                    ci2 = di2 + di5
                    cr3 = dr3 + dr4
                    ci4 = dr4 - dr3
                    cr4 = di3 - di4
                    ci3 = di3 + di4
                    ch(i-1, 1, k) = cc(i-1, k, 1) + cr2 + cr3
                    ch(i, 1, k) = cc(i, k, 1) + ci2 + ci3
                    tr2 = cc(i-1, k, 1) + tr11*cr2 + tr12*cr3
                    ti2 = cc(i, k, 1) + tr11*ci2 + tr12*ci3
                    tr3 = cc(i-1, k, 1) + tr12*cr2 + tr11*cr3
                    ti3 = cc(i, k, 1) + tr12*ci2 + tr11*ci3
                    tr5 = ti11*cr5 + ti12*cr4
                    ti5 = ti11*ci5 + ti12*ci4
                    tr4 = ti12*cr5 - ti11*cr4
                    ti4 = ti12*ci5 - ti11*ci4
                    ch(i-1, 3, k) = tr2 + tr5
                    ch(ic-1, 2, k) = tr2 - tr5
                    ch(i, 3, k) = ti2 + ti5
                    ch(ic, 2, k) = ti5 - ti2
                    ch(i-1, 5, k) = tr3 + tr4
                    ch(ic-1, 4, k) = tr3 - tr4
                    ch(i, 5, k) = ti3 + ti4
                    ch(ic, 4, k) = ti4 - ti3
                end do
            end do
        end if

    end subroutine real_forward_pass_5

    subroutine real_forward_pass_n(ido, iip, l1, idl1, cc, c1, c2, ch, ch2, wa)

        ! Dummy arguments

        integer(ip), intent(in)      :: ido
        integer(ip), intent(in)      :: iip
        integer(ip), intent(in)      :: l1
        integer(ip), intent(in)      :: idl1
        real(wp), intent(out)    :: cc(ido, iip, l1)
        real(wp), intent(inout)  :: c1(ido, l1, iip)
        real(wp), intent(inout)  :: c2(idl1, iip)
        real(wp), intent(inout)  :: ch(ido, l1, iip)
        real(wp), intent(inout)  :: ch2(idl1, iip)
        real(wp), intent(in)     :: wa(*)

        ! Local variables

        integer(ip)         :: iipph, iipp2, idp2, nbd, j
        integer(ip)         :: k, is, idij, i, jc, l, lc
        real(wp)            :: arg, dcp, dsp, ar1, ai1
        real(wp)            :: ar1h, dc2, ds2, ar2, ai2, ar2h


        arg = TWO_PI/iip
        dcp = cos(arg)
        dsp = sin(arg)
        iipph = (iip + 1)/2
        iipp2 = iip + 2
        idp2 = ido + 2
        nbd =(ido - 1)/2

        block_construct: block
            if (ido /= 1) then
                ch2(:, 1) = c2(:, 1)
                ch(1, :, 2:iip) = c1(1, :, 2:iip)
                if (nbd <= l1) then
                    is = -ido
                    do j = 2, iip
                        is = is + ido
                        idij = is
                        do i = 3, ido, 2
                            idij = idij + 2
                            ch(i-1, :, j)=wa(idij-1)*c1(i-1, :, j)+wa(idij)*c1(i, :, j)
                            ch(i, :, j)=wa(idij-1)*c1(i, :, j)-wa(idij)*c1(i-1, :, j)
                        end do
                    end do
                else
                    is = -ido
                    do j = 2, iip
                        is = is + ido
                        do k = 1, l1
                            idij = is
                            ch(2:ido-1:2, k, j) = wa(idij+1:ido-2+idij:2)*c1(2:ido-1 &
                                :2, k, j) + wa(idij+2:ido-1+idij:2)*c1(3:ido:2, k, j)
                            ch(3:ido:2, k, j) = wa(idij+1:ido-2+idij:2)*c1(3:ido:2, k &
                                , j) - wa(idij+2:ido-1+idij:2)*c1(2:ido-1:2, k, j)
                        end do
                    end do
                end if
                if (l1 <= nbd) then
                    do j = 2, iipph
                        jc = iipp2 - j
                        c1(2:ido-1:2, :, j)=ch(2:ido-1:2, :, j)+ch(2:ido-1:2, :, jc)
                        c1(2:ido-1:2, :, jc) = ch(3:ido:2, :, j) - ch(3:ido:2, :, jc)
                        c1(3:ido:2, :, j) = ch(3:ido:2, :, j) + ch(3:ido:2, :, jc)
                        c1(3:ido:2, :, jc) = ch(2:ido-1:2, :, jc) - ch(2:ido-1:2, :, j)
                    end do
                    exit block_construct
                end if
                do j = 2, iipph
                    jc = iipp2 - j
                    c1(2:ido-1:2, :, j) = ch(2:ido-1:2, :, j) + ch(2:ido-1:2, :, jc)
                    c1(2:ido-1:2, :, jc) = ch(3:ido:2, :, j) - ch(3:ido:2, :, jc)
                    c1(3:ido:2, :, j) = ch(3:ido:2, :, j) + ch(3:ido:2, :, jc)
                    c1(3:ido:2, :, jc) = ch(2:ido-1:2, :, jc) - ch(2:ido-1:2, :, j)
                end do
                exit block_construct
            end if
            c2(:, 1) = ch2(:, 1)
        end block block_construct

        do j = 2, iipph
            jc = iipp2 - j
            c1(1, :, j) = ch(1, :, j) + ch(1, :, jc)
            c1(1, :, jc) = ch(1, :, jc) - ch(1, :, j)
        end do

        ar1 = ONE
        ai1 = ZERO
        do l = 2, iipph
            lc = iipp2 - l
            ar1h = dcp*ar1 - dsp*ai1
            ai1 = dcp*ai1 + dsp*ar1
            ar1 = ar1h
            ch2(:, l) = c2(:, 1) + ar1*c2(:, 2)
            ch2(:, lc) = ai1*c2(:, iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j = 3, iipph
                jc = iipp2 - j
                ar2h = dc2*ar2 - ds2*ai2
                ai2 = dc2*ai2 + ds2*ar2
                ar2 = ar2h
                ch2(:, l) = ch2(:, l) + ar2*c2(:, j)
                ch2(:, lc) = ch2(:, lc) + ai2*c2(:, jc)
            end do
        end do

        do j = 2, iipph
            ch2(:, 1) = ch2(:, 1) + c2(:, j)
        end do

        cc(:, 1, :) = ch(:, :, 1)
        cc(ido, 2:(iipph-1)*2:2, :) = transpose(ch(1, :, 2:iipph))
        cc(1, 3:iipph*2-1:2, :) = transpose(ch(1, :, iipp2-2:iipp2-iipph:(-1)))

        if(ido /= 1) then
            if (l1 <= nbd) then
                cc(2:ido-1:2, 3:iipph*2-1:2, :) = &
                    reshape( &
                    source = ch(2:ido-1:2, :, 2:iipph)+ch(2:ido-1:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(idp2-4:idp2-1-ido:(-2), 2:(iipph-1)*2:2, :) = &
                    reshape( &
                    source = ch(2:ido-1:2, :, 2:iipph)-ch(2:ido-1:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(3:ido:2, 3:iipph*2-1:2, :) = &
                    reshape(&
                    source = ch(3:ido:2, :, 2:iipph)+ch(3:ido:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(idp2-3:idp2-ido:(-2), 2:(iipph-1)*2:2, :) = &
                    reshape(&
                    source = ch(3:ido:2, :, iipp2-2:iipp2-iipph:(-1))-ch(3:ido:2, :, 2:iipph), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
            else
                cc(2:ido-1:2, 3:iipph*2-1:2, :) = &
                    reshape( &
                    source = ch(2:ido-1:2, :, 2:iipph)+ch(2:ido-1:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(idp2-4:idp2-1-ido:(-2), 2:(iipph-1)*2:2, :) = &
                    reshape( &
                    source = ch(2:ido-1:2, :, 2:iipph)-ch(2:ido-1:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(3:ido:2, 3:iipph*2-1:2, :) = &
                    reshape( &
                    source = ch(3:ido:2, :, 2:iipph)+ch(3:ido:2, :, iipp2-2:iipp2-iipph:(-1)), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
                cc(idp2-3:idp2-ido:(-2), 2:(iipph-1)*2:2, :) = &
                    reshape( &
                    source = ch(3:ido:2, :, iipp2-2:iipp2-iipph:(-1))-ch(3:ido:2, :, 2:iipph), &
                    shape = [(ido-1)/2, iipph-1, l1], &
                    order = [1, 3, 2] &
                    )
            end if
        end if

    end subroutine real_forward_pass_n

end module type_FastFourierTransform
