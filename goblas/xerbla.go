package goblas

// \brief \b Xerbla
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Xerbla( SRNAME, INFO)
//
//       .. Scalar Arguments ..
//       CHARACTER//(//)      SRNAME
//       INTEGER            INFO
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Xerbla  is an error handler for the LAPACK routines.
// It is called by an LAPACK routine if an input parameter has an
// invalid value.  A message is printed and execution stops.
//
// Installers may consider modifying the STOP statement in order to
// call system-specific exception-handling facilities.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SRNAME
// \verbatim
//          SRNAME is CHARACTER//(//)
//          The name of the routine which called Xerbla.
// \endverbatim
//
// \param[in] INFO
// \verbatim
//          INFO is INTEGER
//          The position of the invalid parameter in the parameter list
//          of the calling routine.
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
// \ingroup aux_blas
//
//  =====================================================================
func Xerbla(srname *[]byte, info *int) {
	//*
	//*  -- Reference BLAS level1 routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     December 2016
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*
	//* =====================================================================
	//*
	//*     .. Intrinsic Functions ..
	//*     ..
	//*     .. Executable Statements ..
	//*
	WRITE(6, *func() *[]byte {
		y := []byte(" ** on entry to %s parameter number %2d had an illegal value\n")
		return &y
	}(), (*srname)[1:LEN_TRIM(srname)-1], (*info))
	//*
	panic("")
	//*
	//*
	//*     End of Xerbla
	//*
}
