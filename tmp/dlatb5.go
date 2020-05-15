package goblas

import 

// Dlatb5 sets parameters for the matrix generator based on the type
// of matrix to be generated.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlatb5( path, imat, n, _type, kl, ku, anorm, mode,
//                          cndnum, dist)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   anorm, cndnum
//       inTEGER            imat, kl, ku, mode, N
//       CHARACTER          dist, _type
//       CHARACTER*3        path
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlatb5 sets parameters for the matrix generator based on the type
// of matrix to be generated.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The lapACK path name.
// \endverbatim
//
// \param[in] imat
// \verbatim
//          imat is inTEGER
//          An integer key describing which matrix to generate for this
//          path.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of rows and columns in the matrix to be generated.
// \endverbatim
//
// \param[out] _type
// \verbatim
//          _type is CHARACTER*1
//          The type of the matrix to be generated:
//          = 'S':  symmetric matrix
//          = 'P':  symmetric positive (semi)definite matrix
//          = 'N':  nonsymmetric matrix
// \endverbatim
//
// \param[out] kl
// \verbatim
//          kl is inTEGER
//          The lower band width of the matrix to be generated.
// \endverbatim
//
// \param[out] ku
// \verbatim
//          ku is inTEGER
//          The upper band width of the matrix to be generated.
// \endverbatim
//
// \param[out] anorm
// \verbatim
//          anorm is DOUBLE PRECISION
//          The desired norm of the matrix to be generated.  The diagonal
//          matrix of singular values or eigenvalues is scaled by this
//          value.
// \endverbatim
//
// \param[out] mode
// \verbatim
//          mode is inTEGER
//          A key indicating how to choose the vector of eigenvalues.
// \endverbatim
//
// \param[out] cndnum
// \verbatim
//          cndnum is DOUBLE PRECISION
//          The desired condition number.
// \endverbatim
//
// \param[out] dist
// \verbatim
//          dist is CHARACTER*1
//          The type of distribution to be used by the random number
//          generator.
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
// \date December 2016
//
// \ingroup double_lin
//
//  =====================================================================
func Dlatb5(path *[]byte, imat *int, n *int, _type *byte, kl *int, ku *int, anorm *float64, mode *int, cndnum *float64, dist *byte) {
	shrink := new(float64)
	tenth := new(float64)
	one := new(float64)
	two := new(float64)
	badc1 := new(float64)
	badc2 := new(float64)
	eps := new(float64)
	large := new(float64)
	small := new(float64)
	first := new(bool)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*shrink) = 0.25
	(*tenth) = 0.1e+0
	(*one) = 1.0e+0
	(*two) = 2.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Save statement ..
	//     ..
	//     .. Data statements ..
	(*first) = true
	//     ..
	//     .. Executable Statements ..
	//
	//     Set some _constants for use in the subroutine.
	//
	if *first {
		(*first) = false
		(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("Precision"); return &y}()))
		(*badc2) = (*tenth) / (*eps)
		(*badc1) = (SQRt((*badc2)))
		(*small) = (*Dlamch(func() *[]byte {y :=[]byte("Safe minimum"); return &y}()))
		(*large) = (*one) / (*small)
		//
		//        If it looks like we're on a Cray, take the square root of
		//        small and large to avoid overflow and underflow problems.
		//
		Dlabad(small, large)
		(*small) = (*shrink) * ((*small) / (*eps))
		(*large) = (*one) / (*small)
	}
	//
	(*c2) = (*(path))[1]
	//
	//     Set some parameters
	//
	(*(dist)) = 'S'
	(*(mode)) = 3
	//
	//     Set _type, the type of matrix to be generated.
	//
	(*(_type)) = (*c2)[0]
	//
	//     Set the lower and upper bandwidths.
	//
	if (*(imat)) == 1 {
		(*(kl)) = 0
	} else {
		(*(kl)) = (MAX((*(n))-1, 0))
	}
	(*(ku)) = (*(kl))
	//
	//     Set the condition number and norm.etc
	//
	if (*(imat)) == 3 {
		(*(cndnum)) = 1.0e12
		(*(mode)) = 2
	} else if (*(imat)) == 4 {
		(*(cndnum)) = 1.0e12
		(*(mode)) = 1
	} else if (*(imat)) == 5 {
		(*(cndnum)) = 1.0e12
		(*(mode)) = 3
	} else if (*(imat)) == 6 {
		(*(cndnum)) = (*badc1)
	} else if (*(imat)) == 7 {
		(*(cndnum)) = (*badc2)
	} else {
		(*(cndnum)) = (*two)
	}
	//
	if (*(imat)) == 8 {
		(*(anorm)) = (*small)
	} else if (*(imat)) == 9 {
		(*(anorm)) = (*large)
	} else {
		(*(anorm)) = (*one)
	}
	//
	if (*(n)) <= 1 {
		(*(cndnum)) = (*one)
	}
	//
	return
	//
	//     End of Dlatb5
	//
}
