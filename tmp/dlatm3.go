package goblas

//    dlatm3 returns the (isub,jsub) entry of a random matrix of
//    dimension (m, N) described by the other parameters. (isub,jsub)
//    is the final position of the (I,j) entry after pivoting
//    according to ipvtng and iwork. dlatm3 is called by the
//    DLAtmR routine in order to build random test matrices. No error
//    checking on parameters is done, because this routine is called in
//    a tight loop by DLAtmR which has already checked the parameters.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION dlatm3( m, n, i, j, isub, jsub, kl, ku,
//                        idist, iseed, d, igrade, DL, DR, ipvtng, iwork,
//                        sparse)
//
//       .. Scalar Arguments ..
//
//       inTEGER            i, idist, igrade, ipvtng, isub, j, jsub, kl,
//      $                   ku, m, N
//       DOUBLE PRECISION   sparse
//       ..
//
//       .. Array Arguments ..
//
//       inTEGER            iseed( 4), iwork(*)
//       DOUBLE PRECISION   d(*), DL(*), dr(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    dlatm3 returns the (isub,jsub) entry of a random matrix of
//    dimension (m, N) described by the other parameters. (isub,jsub)
//    is the final position of the (I,j) entry after pivoting
//    according to ipvtng and iwork. dlatm3 is called by the
//    DLAtmR routine in order to build random test matrices. No error
//    checking on parameters is done, because this routine is called in
//    a tight loop by DLAtmR which has already checked the parameters.
//
//    Use of dlatm3 differs from SLATm2 in the order in which the random
//    number generator is called to fill in random matrix entries.
//    With Dlatm2, the generator is called to fill in the pivoted matrix
//    columnwise. With dlatm3, the generator is called to fill in the
//    matrix columnwise, after which it is pivoted. Thus, dlatm3 can
//    be used to _construct random matrices which differ only in their
//    order of rows and/or columns. Dlatm2 is used to _construct band
//    matrices while avoiding calling the random number generator for
//    entries outside the band (and therefore generating random numbers
//    in different orders for different pivot orders).
//
//    The matrix whose (isub,jsub) entry is returned is _constructed as
//    follows (this routine only computes one entry):
//
//      If isub is outside (1..M) or jsub is outside (1..N), return zero
//         (this is convenient for generating matrices in band format).
//
//      Generate a matrix A with random entries of distribution idist.
//
//      Set the diagonal to D.
//
//      Grade the matrix, if desired, from the left (by DL) and/or
//         from the right (by DR or DL) as specified by igrade.
//
//      Permute, if desired, the rows and/or columns as specified by
//         ipvtng and iwork.
//
//      Band the matrix to have lower bandwidth kl and upper
//         bandwidth ku.
//
//      Set random entries to zero as specified by sparse.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//           Number of rows of matrix. Not modified.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//           Number of columns of matrix. Not modified.
// \endverbatim
//
// \param[in] I
// \verbatim
//          I is inTEGER
//           Row of unpivoted entry to be returned. Not modified.
// \endverbatim
//
// \param[in] J
// \verbatim
//          J is inTEGER
//           Column of unpivoted entry to be returned. Not modified.
// \endverbatim
//
// \param[in,out] isub
// \verbatim
//          isub is inTEGER
//           Row of pivoted entry to be returned. Changed on exit.
// \endverbatim
//
// \param[in,out] jsub
// \verbatim
//          jsub is inTEGER
//           Column of pivoted entry to be returned. Changed on exit.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is inTEGER
//           Lower bandwidth. Not modified.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is inTEGER
//           Upper bandwidth. Not modified.
// \endverbatim
//
// \param[in] idist
// \verbatim
//          idist is inTEGER
//           On entry, idist specifies the type of distribution to be
//           used to generate a random matrix .
//           1 => UNIFORM( 0, 1)
//           2 => UNIFORM( -1, 1)
//           3 => normaL( 0, 1)
//           Not modified.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is inTEGER array of dimension ( 4)
//           Seed for random number generator.
//           Changed on exit.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array of dimension ( Min( I, J))
//           diagonal entries of matrix. Not modified.
// \endverbatim
//
// \param[in] igrade
// \verbatim
//          igrade is inTEGER
//           Specifies grading of matrix as follows:
//           0  => no grading
//           1  => matrix premultiplied by diag( DL)
//           2  => matrix postmultiplied by diag( DR)
//           3  => matrix premultiplied by diag( DL) and
//                         postmultiplied by diag( DR)
//           4  => matrix premultiplied by diag( DL) and
//                         postmultiplied by inv( diag( DL))
//           5  => matrix premultiplied by diag( DL) and
//                         postmultiplied by diag( DL)
//           Not modified.
// \endverbatim
//
// \param[in] DL
// \verbatim
//          DL is DOUBLE PRECISION array ( I or j, as appropriate)
//           Left scale factors for grading matrix.  Not modified.
// \endverbatim
//
// \param[in] DR
// \verbatim
//          DR is DOUBLE PRECISION array ( I or j, as appropriate)
//           Right scale factors for grading matrix.  Not modified.
// \endverbatim
//
// \param[in] ipvtng
// \verbatim
//          ipvtng is inTEGER
//           On entry specifies pivoting permutations as follows:
//           0 => none.
//           1 => row pivoting.
//           2 => column pivoting.
//           3 => full pivoting, i.e., on both sides.
//           Not modified.
// \endverbatim
//
// \param[in] iwork
// \verbatim
//          iwork is inTEGER array ( I or j, as appropriate)
//           This array specifies the permutation used. The
//           row (or column) originally in position K is in
//           position iwork( K) after pivoting.
//           This differs from iwork for Dlatm2. Not modified.
// \endverbatim
//
// \param[in] sparse
// \verbatim
//          sparse is DOUBLE PRECISION between 0. and 1.
//           On entry specifies the sparsity of the matrix
//           if sparse matrix is to be generated.
//           sparse should lie between 0 and 1.
//           A uniform ( 0, 1) random number x is generated and
//           compared to sparse; if x is larger the matrix entry
//           is unchanged and if x is smaller the entry is set
//           to zero. Thus on the average a fraction sparse of the
//           entries will be set to zero.
//           Not modified.
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG Ltd.
//
// \date June 2016
//
// \ingroup double_matgen
//
//  =====================================================================
func dlatm3(m *int, n *int, i *int, j *int, isub *int, jsub *int, kl *int, ku *int, idist *int, iseed *[]int, d *[]float64, igrade *int, dl *[]float64, dr *[]float64, ipvtng *int, iwork *[]int, sparse *float64) (dlatm3Return *float64) {
	dlatm3Return = new(float64)
	zero := new(float64)
	temp := new(float64)
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2016
	//
	//     .. Scalar Arguments ..
	//
	//     ..
	//
	//     .. Array Arguments ..
	//
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	//
	(*zero) = 0.0
	//     ..
	//
	//     .. Local Scalars ..
	//
	//     ..
	//
	//     .. External Functions ..
	//
	//     ..
	//
	//-----------------------------------------------------------------------
	//
	//     .. Executable Statements ..
	//
	//
	//     Check for I and J in range
	//
	if (*(i)) < 1 || (*(i)) > (*(m)) || (*(j)) < 1 || (*(j)) > (*(n)) {
		(*(isub)) = (*(i))
		(*(jsub)) = (*(j))
		(*(dlatm3Return)) = (*zero)
		return
	}
	//
	//     Compute subscripts depending on ipvtng
	//
	if (*(ipvtng)) == 0 {
		(*(isub)) = (*(i))
		(*(jsub)) = (*(j))
	} else if (*(ipvtng)) == 1 {
		(*(isub)) = (*(iwork))[(*(i))-1]
		(*(jsub)) = (*(j))
	} else if (*(ipvtng)) == 2 {
		(*(isub)) = (*(i))
		(*(jsub)) = (*(iwork))[(*(j))-1]
	} else if (*(ipvtng)) == 3 {
		(*(isub)) = (*(iwork))[(*(i))-1]
		(*(jsub)) = (*(iwork))[(*(j))-1]
	}
	//
	//     Check for banding
	//
	if (*(jsub)) > (*(isub))+(*(ku)) || (*(jsub)) < (*(isub))-(*(kl)) {
		(*(dlatm3Return)) = (*zero)
		return
	}
	//
	//     Check for sparsity
	//
	if (*(sparse)) > (*zero) {
		if (*Dlaran((iseed))) < (*(sparse)) {
			(*(dlatm3Return)) = (*zero)
			return
		}
	}
	//
	//     Compute entry and grade it according to igrade
	//
	if (*(i)) == (*(j)) {
		(*temp) = (*(d))[(*(i))-1]
	} else {
		(*temp) = (*Dlarnd((idist), (iseed)))
	}
	if (*(igrade)) == 1 {
		(*temp) = (*temp) * (*(dl))[(*(i))-1]
	} else if (*(igrade)) == 2 {
		(*temp) = (*temp) * (*(dr))[(*(j))-1]
	} else if (*(igrade)) == 3 {
		(*temp) = (*temp) * (*(dl))[(*(i))-1] * (*(dr))[(*(j))-1]
	} else if (*(igrade)) == 4 && (*(i)) != (*(j)) {
		(*temp) = (*temp) * (*(dl))[(*(i))-1] / (*(dl))[(*(j))-1]
	} else if (*(igrade)) == 5 {
		(*temp) = (*temp) * (*(dl))[(*(i))-1] * (*(dl))[(*(j))-1]
	}
	(*(dlatm3Return)) = (*temp)
	return
	//
	//     End of dlatm3
	//
}
