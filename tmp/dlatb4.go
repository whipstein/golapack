package goblas

import 

// Dlatb4 sets parameters for the matrix generator based on the type of
// matrix to be generated.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlatb4( path, imat, m, n, _type, kl, ku, anorm, mode,
//                          cndnum, dist)
//
//       .. Scalar Arguments ..
//       CHARACTER          dist, _type
//       CHARACTER*3        path
//       intEGER            imat, kl, ku, m, mode, N
//       DOUBLE PRECISION   anorm, cndnum
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlatb4 sets parameters for the matrix generator based on the type of
// matrix to be generated.
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
//          imat is intEGER
//          An integer key describing which matrix to generate for this
//          path.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows in the matrix to be generated.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns in the matrix to be generated.
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
//          kl is intEGER
//          The lower band width of the matrix to be generated.
// \endverbatim
//
// \param[out] ku
// \verbatim
//          ku is intEGER
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
//          mode is intEGER
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
func Dlatb4(path *[]byte, imat *int, m *int, n *int, _type *byte, kl *int, ku *int, anorm *float64, mode *int, cndnum *float64, dist *byte) {
	SHRInk := new(float64)
	tenth := new(float64)
	one := new(float64)
	two := new(float64)
	first := new(bool)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	mat := new(int)
	badc1 := new(float64)
	badc2 := new(float64)
	eps := new(float64)
	large := new(float64)
	small := new(float64)
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
		(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("Precision"); return &y }()))
		(*badc2) = (*tenth) / (*eps)
		(*badc1) = (SQRt((*badc2)))
		(*small) = (*Dlamch(func() *[]byte {y :=[]byte("Safe minimum"); return &y }()))
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
	//     Set some parameters we don't plan to change.
	//
	(*(dist)) = 'S'
	(*(mode)) = 3
	//
	if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("QR"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("LQ"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("QL"); return &y }()) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("RQ"); return &y }())) {
		//
		//        xQR, xLq, xQL, xRQ:  Set parameters to generate a general
		//                             M x N matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the lower and upper bandwidths.
		//
		if (*(imat)) == 1 {
			(*(kl)) = 0
			(*(ku)) = 0
		} else if (*(imat)) == 2 {
			(*(kl)) = 0
			(*(ku)) = (MAX((*(n))-1, 0))
		} else if (*(imat)) == 3 {
			(*(kl)) = (MAX((*(m))-1, 0))
			(*(ku)) = 0
		} else {
			(*(kl)) = (MAX((*(m))-1, 0))
			(*(ku)) = (MAX((*(n))-1, 0))
		}
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 5 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 6 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 7 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 8 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("GE"); return &y }()) {
		//
		//        xGE:  Set parameters to generate a general M x N matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the lower and upper bandwidths.
		//
		if (*(imat)) == 1 {
			(*(kl)) = 0
			(*(ku)) = 0
		} else if (*(imat)) == 2 {
			(*(kl)) = 0
			(*(ku)) = (MAX((*(n))-1, 0))
		} else if (*(imat)) == 3 {
			(*(kl)) = (MAX((*(m))-1, 0))
			(*(ku)) = 0
		} else {
			(*(kl)) = (MAX((*(m))-1, 0))
			(*(ku)) = (MAX((*(n))-1, 0))
		}
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 8 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 9 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 10 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 11 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("GB"); return &y }()) {
		//
		//        xGB:  Set parameters to generate a general banded matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 5 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 6 {
			(*(cndnum)) = (*tenth) * (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 7 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 8 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("GT"); return &y }()) {
		//
		//        xGT:  Set parameters to generate a general tridiagonal matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the lower and upper bandwidths.
		//
		if (*(imat)) == 1 {
			(*(kl)) = 0
		} else {
			(*(kl)) = 1
		}
		(*(ku)) = (*(kl))
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 3 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 4 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 5 || (*(imat)) == 11 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 6 || (*(imat)) == 12 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("PO"); return &y }())) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("PP"); return &y }())) {
		//
		//        xPO, xPP: Set parameters to generate a
		//        symmetric positive definite matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = (*c2)[0]
		//
		//        Set the lower and upper bandwidths.
		//
		if (*(imat)) == 1 {
			(*(kl)) = 0
		} else {
			(*(kl)) = (MAX((*(n))-1, 0))
		}
		(*(ku)) = (*(kl))
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 6 {
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
		//
	} else if (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("SY"); return &y }())) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("SP"); return &y }())) {
		//
		//        xSY, xSP: Set parameters to generate a
		//        symmetric matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = (*c2)[0]
		//
		//        Set the lower and upper bandwidths.
		//
		if (*(imat)) == 1 {
			(*(kl)) = 0
		} else {
			(*(kl)) = (MAX((*(n))-1, 0))
		}
		(*(ku)) = (*(kl))
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 7 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 8 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 9 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 10 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("PB"); return &y }()) {
		//
		//        xPB:  Set parameters to generate a symmetric band matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'P'
		//
		//        Set the norm and condition number.
		//
		if (*(imat)) == 5 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 6 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 7 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 8 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("PT"); return &y }()) {
		//
		//        xPT:  Set parameters to generate a symmetric positive definite
		//        tridiagonal matrix.
		//
		(*(_type)) = 'P'
		if (*(imat)) == 1 {
			(*(kl)) = 0
		} else {
			(*(kl)) = 1
		}
		(*(ku)) = (*(kl))
		//
		//        Set the condition number and norm.
		//
		if (*(imat)) == 3 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 4 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 5 || (*(imat)) == 11 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 6 || (*(imat)) == 12 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TR"); return &y }())) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TP"); return &y }())) {
		//
		//        xTR, xTP:  Set parameters to generate a triangular matrix
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the lower and upper bandwidths.
		//
		(*mat) = (ABS((*(imat))))
		if (*mat) == 1 || (*mat) == 7 {
			(*(kl)) = 0
			(*(ku)) = 0
		} else if (*(imat)) < 0 {
			(*(kl)) = (MAX((*(n))-1, 0))
			(*(ku)) = 0
		} else {
			(*(kl)) = 0
			(*(ku)) = (MAX((*(n))-1, 0))
		}
		//
		//        Set the condition number and norm.
		//
		if (*mat) == 3 || (*mat) == 9 {
			(*(cndnum)) = (*badc1)
		} else if (*mat) == 4 {
			(*(cndnum)) = (*badc2)
		} else if (*mat) == 10 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*mat) == 5 {
			(*(anorm)) = (*small)
		} else if (*mat) == 6 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TB"); return &y }()) {
		//
		//        xTB:  Set parameters to generate a triangular band matrix.
		//
		//        Set _type, the type of matrix to be generated.
		//
		(*(_type)) = 'N'
		//
		//        Set the norm and condition number.
		//
		if (*(imat)) == 2 || (*(imat)) == 8 {
			(*(cndnum)) = (*badc1)
		} else if (*(imat)) == 3 || (*(imat)) == 9 {
			(*(cndnum)) = (*badc2)
		} else {
			(*(cndnum)) = (*two)
		}
		//
		if (*(imat)) == 4 {
			(*(anorm)) = (*small)
		} else if (*(imat)) == 5 {
			(*(anorm)) = (*large)
		} else {
			(*(anorm)) = (*one)
		}
	}
	if (*(n)) <= 1 {
		(*(cndnum)) = (*one)
	}
	//
	return
	//
	//     End of Dlatb4
	//
}
