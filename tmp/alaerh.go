package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Alaerh is an error handler for the lapACK routines.  It prints the
// header if this is the first error message and prints the error code
// and form of recovery, if any.  The character evaluations in this
// routine may make it slow, but it should not be called once the lapACK
// routines are fully debugged.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Alaerh( path, subnam, info, infoe, opts, m, n, kl, ku,
//                          n5, imat, nfail, nerrs, nout)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        path
//       CHARACTER*(*)    subnam
//       CHARACTER*(*)    opts
//       intEGER            imat, info, infoe, kl, ku, m, n, n5, nerrs,
//      $                   nfail, nout
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Alaerh is an error handler for the lapACK routines.  It prints the
// header if this is the first error message and prints the error code
// and form of recovery, if any.  The character evaluations in this
// routine may make it slow, but it should not be called once the lapACK
// routines are fully debugged.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The lapACK path name of subroutine subnam.
// \endverbatim
//
// \param[in] subnam
// \verbatim
//          subnam is CHARACTER*(*)
//          The name of the subroutine that returned an error code.
// \endverbatim
//
// \param[in] info
// \verbatim
//          info is intEGER
//          The error code returned from routine subnam.
// \endverbatim
//
// \param[in] infoe
// \verbatim
//          infoe is intEGER
//          The expected error code from routine subnam, if subnam were
//          error-free.  If infoe = 0, an error message is printed, but
//          if infoe.NE.0, we assume only the return code info is wrong.
// \endverbatim
//
// \param[in] opts
// \verbatim
//          opts is CHARACTER*(*)
//          The character options to the subroutine subnam, concatenated
//          into a single character string.  For example, uplo = 'U',
//          trans = 'T', and diag = 'N' for a triangular routine would
//          be specified as opts = 'UTN'.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The matrix row dimension.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The matrix column dimension.  Accessed only if path = xGE or
//          xGB.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is intEGER
//          The number of sub-diagonals of the matrix.  Accessed only if
//          path = xGB, xPB, or xTB.  Also used for nrhs for path = xLS.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//          The number of super-diagonals of the matrix.  Accessed only
//          if path = xGB.
// \endverbatim
//
// \param[in] n5
// \verbatim
//          n5 is intEGER
//          A fifth integer parameter, may be the blocksize nb or the
//          number of right hand sides nrhs.
// \endverbatim
//
// \param[in] imat
// \verbatim
//          imat is intEGER
//          The matrix type.
// \endverbatim
//
// \param[in] nfail
// \verbatim
//          nfail is intEGER
//          The number of prior tests that did not pass the threshold;
//          used to determine if the header should be printed.
// \endverbatim
//
// \param[in,out] nerrs
// \verbatim
//          nerrs is intEGER
//          On entry, the number of errors already detected; used to
//          determine if the header should be printed.
//          On exit, nerrs is increased by 1.
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number on which results are to be printed.
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
// \ingroup aux_lin
//
//  =====================================================================
func Alaerh(path *[]byte, subnam *[]byte, info *int, infoe *int, opts *[]byte, m *int, n *int, kl *int, ku *int, n5 *int, imat *int, nfail *int, nerrs *int, nout *int) {
	uplo := new(byte)
	p2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	c3 := func() *[]byte {
		arr := make([]byte, 3)
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
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*(info)) == 0 {
		return
	}
	(*p2)[0] = (*(path))[1]
	(*c3)[0] = (*(subnam))[3]
	//
	//     Print the header if this is the first error message.
	//
	if (*(nfail)) == 0 && (*(nerrs)) == 0 {
		if (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }())) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }())) {
			Aladhd((nout), (path))
		} else {
			Alahd((nout), (path))
		}
	}
	(*(nerrs)) = (*(nerrs)) + 1
	//
	//     Print the message detailing the error and form of recovery,
	//     if any.
	//
	if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GE"); return &y }()) {
		//
		//        xGE:  General matrices
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> M =%5d, N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(m)), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s=%5d for M=%5d, N=%5d, nb=%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d for N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', trans='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', trans='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s=%5d for N=%5d, nb=%4d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(n5)), (*(imat)))
			//
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for M =%5d, N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(imat)))
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for NORM = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(imat)))
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("LS "); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> trans = '%c', M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
			//
		} else if (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("LSX"); return &y }())) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("LSS"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> trans = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GB"); return &y }()) {
		//
		//        xGB:  General band matrices
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> M = %5d, N =%5d, kl =%5d, ku =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(m)), (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> M = %5d, N =%5d, kl =%5d, ku =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> N =%5d, kl =%5d, ku =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> N =%5d, kl =%5d, ku =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', trans='%c', N=%5d, kl=%5d, ku=%5d, nrhs=%4d, type %1d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', trans='%c', N=%5d, kl=%5d, ku=%5d, nrhs=%4d, type %1d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> M = %5d, N =%5d, kl =%5d, ku =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(ku)), (*(imat)))
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> NORM ='%c', N =%5d, kl =%5d, ku =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(kl)), (*(ku)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> trans='%c', N =%5d, kl =%5d, ku =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(kl)), (*(ku)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GT"); return &y }()) {
		//
		//        xGT:  General tridiagonal matrices
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d for N=%5d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d for N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', trans='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', trans='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for NORM = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> trans = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PO"); return &y }()) {
		//
		//        xPO:  Symmetric or Hermitian positive definite matrices
		//
		(*uplo) = (*(opts))[0]
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			//
		} else if (Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }())) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for uplo = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PS"); return &y }()) {
		//
		//        xPS:  Symmetric or Hermitian positive semi-definite matrices
		//
		(*uplo) = (*(opts))[0]
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*(infoe)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*(infoe)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam)), (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam)), (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			//
		} else if (Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmT"); return &y }())) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for uplo = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam)), (*(info)), (*uplo), (*(m)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam)), (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SY"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SR"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SK"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HE"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HR"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HK"); return &y }()) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HA"); return &y }())) {
		//
		//        xSY: symmetric indefinite matrices
		//             with partial (Bunch-Kaufman) pivoting;
		//        xSR: symmetric indefinite matrices
		//             with rook (bounded Bunch-Kaufman) pivoting;
		//        xSK: symmetric indefinite matrices
		//             with rook (bounded Bunch-Kaufman) pivoting,
		//             new storage format;
		//        xHE: Hermitian indefinite matrices
		//             with partial (Bunch-Kaufman) pivoting.
		//        xHR: Hermitian indefinite matrices
		//             with rook (bounded Bunch-Kaufman) pivoting;
		//        xHK: Hermitian indefinite matrices
		//             with rook (bounded Bunch-Kaufman) pivoting,
		//             new storage format;
		//        xHA: Hermitian matrices
		//             Aasen Algorithm
		//
		(*uplo) = (*(opts))[0]
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 2; return &y }(), c3, func() *[]byte {y := []byte("SV"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) || Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for uplo = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PP"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SP"); return &y }()) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HP"); return &y }())) {
		//
		//        xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices
		//
		(*uplo) = (*(opts))[0]
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(m)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d for uplo = '%c', N =%5d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', uplo='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) || Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d for uplo = '%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PB"); return &y }()) {
		//
		//        xPB:  Symmetric (Hermitian) positive definite band matrix
		//
		(*uplo) = (*(opts))[0]
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo = '%c', N =%5d, kd =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(m)), (*(kl)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, kd =%5d, nb =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(kl)), (*(n5)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> uplo='%c', N =%5d, kd =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*uplo), (*(n)), (*(kl)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s=%5d\n ==> uplo = '%c', N =%5d, kd =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(n)), (*(kl)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', uplo='%c', N=%5d, kd=%5d, nrhs=%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(opts))[1], (*(n)), (*(kl)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d\n ==> fact='%c', uplo='%c', N=%5d, kd=%5d, nrhs=%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(n)), (*(kl)), (*(n5)), (*(imat)))
			}
			//
		} else if (Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }())) || (Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }())) {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo = '%c', N =%5d, kd =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(kl)), (*(imat)))
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> uplo = '%c', N =%5d, kd =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*uplo), (*(m)), (*(kl)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PT"); return &y }()) {
		//
		//        xPT:  Positive definite tridiagonal matrices
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("TRF"); return &y }()) {
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d for N=%5d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(imat)))
			}
			if (*(info)) != 0 {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" ==> Doing only the condition estimate for this case\n"); return &y }())
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SV "); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d for N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("SVX"); return &y }()) {
			//
			if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> fact='%c', N =%5d, nrhs =%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(opts))[0], (*(n)), (*(n5)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s=%5d, fact='%c', N=%5d, nrhs=%4d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(n)), (*(n5)), (*(imat)))
			}
			//
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			//
			if (*blas.Lsame(&((*(subnam))[0]), func() *byte {y := byte('S'); return &y }())) || (*blas.Lsame(&((*(subnam))[0]), func() *byte {y := byte('D'); return &y }())) {
				WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(imat)))
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y := []byte(" *** Error code from %s =%5d for NORM = '%c', N =%5d, type %2d\n")
					return &y
				}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(imat)))
			}
			//
		} else {
			//
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> trans = '%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TR"); return &y }()) {
		//
		//        xTR:  Triangular matrix
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', diag ='%c', N =%5d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(m)), (*(n5)), (*(imat)))
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> NORM='%c', uplo ='%c', diag='%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LATRS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', NOrmin='%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(opts))[3], (*(m)), (*(imat)))
		} else {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TP"); return &y }()) {
		//
		//        xTP:  Triangular packed matrix
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("tri"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', diag ='%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(m)), (*(imat)))
		} else if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> NORM='%c', uplo ='%c', diag='%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LATPS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', NOrmin='%c', N =%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(opts))[3], (*(m)), (*(imat)))
		} else {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', N =%5d, nrhs =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TB"); return &y }()) {
		//
		//        xTB:  Triangular band matrix
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("CON"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> NORM='%c', uplo ='%c', diag='%c', N=%5d, kd=%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(kl)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LATBS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', NOrmin='%c', N=%5d, kd=%5d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(opts))[3], (*(m)), (*(kl)), (*(imat)))
		} else {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s =%5d\n ==> uplo='%c', trans='%c', diag='%c', N=%5d, kd=%5d, nrhs=%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(opts))[0], (*(opts))[1], (*(opts))[2], (*(m)), (*(kl)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QR"); return &y }()) {
		//
		//        xQR:  QR factorization
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("qrs"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for M =%5d, N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("LQ"); return &y }()) {
		//
		//        xLQ:  LQ factorization
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("LQS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for M =%5d, N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QL"); return &y }()) {
		//
		//        xQL:  QL factorization
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("QLS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for M =%5d, N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("RQ"); return &y }()) {
		//
		//        xRQ:  RQ factorization
		//
		if Lsamen(func() *int {y := 3; return &y }(), c3, func() *[]byte {y := []byte("RQS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d\n ==> M =%5d, N =%5d, nrhs =%4d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(kl)), (*(n5)), (*(imat)))
		} else if Lsamen(func() *int {y := 5; return &y }(), &((*(subnam))[1]), func() *[]byte {y := []byte("LAtmS"); return &y }()) {
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d for M =%5d, N =%5d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("LU"); return &y }()) {
		//
		if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> M =%5d, N =%5d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(m)), (*(n)), (*(n5)), (*(imat)))
		} else {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** Error code from %s=%5d for M=%5d, N=%5d, nb=%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n)), (*(n5)), (*(imat)))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("CH"); return &y }()) {
		//
		if (*(info)) != (*(infoe)) && (*(infoe)) != 0 {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" *** %s returned with info =%5d instead of %2d\n ==> N =%5d, nb =%4d, type %2d\n")
				return &y
			}(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(infoe)), (*(m)), (*(n5)), (*(imat)))
		} else {
			WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s=%5d for N=%5d, nb=%4d, type %2d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)), (*(m)), (*(n5)), (*(imat)))
		}
		//
	} else {
		//
		//        Print a generic message if the path is unknown.
		//
		WRITE((*(nout)), *func() *[]byte {y := []byte(" *** Error code from %s =%5d\n"); return &y }(), (*(subnam))[1:lenTrim((subnam))-1], (*(info)))
	}
	//
	//     Description of error message (alphabetical, left to right)
	//
	//     subnam, info, fact, n, nrhs, imat
	//
	//
	//     subnam, info, fact, trans, n, kl, ku, nrhs, imat
	//
	//
	//     subnam, info, fact, trans, n, nrhs, imat
	//
	//
	//     subnam, info, fact, uplo, n, kd, nrhs, imat
	//
	//
	//     subnam, info, fact, uplo, n, nrhs, imat
	//
	//
	//     subnam, info, infoe, fact, n, nrhs, imat
	//
	//
	//     subnam, info, infoe, fact, trans, n, kl, ku, nrhs, imat
	//
	//
	//     subnam, info, infoe, fact, trans, n, nrhs, imat
	//
	//
	//     subnam, info, infoe, fact, uplo, n, kd, nrhs, imat
	//
	//
	//     subnam, info, infoe, fact, uplo, n, nrhs, imat
	//
	//
	//     subnam, info, infoe, m, n, kl, ku, nb, imat
	//
	//
	//     subnam, info, infoe, m, n, nb, imat
	//
	//
	//     subnam, info, infoe, n, imat
	//
	//
	//     subnam, info, infoe, n, kl, ku, nrhs, imat
	//
	//
	//     subnam, info, infoe, n, nb, imat
	//
	//
	//     subnam, info, infoe, n, nrhs, imat
	//
	//
	//     subnam, info, infoe, uplo, n, imat
	//
	//
	//     subnam, info, infoe, uplo, n, kd, nb, imat
	//
	//
	//     subnam, info, infoe, uplo, n, kd, nrhs, imat
	//
	//
	//     subnam, info, infoe, uplo, n, nb, imat
	//
	//
	//     subnam, info, infoe, uplo, n, nrhs, imat
	//
	//
	//     subnam, info, m, n, imat
	//
	//
	//     subnam, info, m, n, kl, ku, imat
	//
	//
	//     subnam, info, m, n, kl, ku, nb, imat
	//
	//
	//     subnam, info, m, n, nb, imat
	//
	//
	//     subnam, info, m, n, nrhs, nb, imat
	//
	//
	//     subnam, info, n, imat
	//
	//
	//     subnam, info, n, kl, ku, nrhs, imat
	//
	//
	//     subnam, info, n, nb, imat
	//
	//
	//     subnam, info, n, nrhs, imat
	//
	//
	//     subnam, info, NORM, n, imat
	//
	//
	//     subnam, info, NORM, n, kl, ku, imat
	//
	//
	//     subnam, info, NORM, uplo, diag, n, imat
	//
	//
	//     subnam, info, NORM, uplo, diag, n, kd, imat
	//
	//
	//     subnam, info, trans, m, n, nrhs, nb, imat
	//
	//
	//     subnam, info, trans, n, kl, ku, nrhs, imat
	//
	//
	//     subnam, info, trans, n, nrhs, imat
	//
	//
	//     subnam, info, uplo, diag, n, imat
	//
	//
	//     subnam, info, uplo, diag, n, nb, imat
	//
	//
	//     subnam, info, uplo, n, imat
	//
	//
	//     subnam, info, uplo, n, kd, imat
	//
	//
	//     subnam, info, uplo, n, kd, nb, imat
	//
	//
	//     subnam, info, uplo, n, kd, nrhs, imat
	//
	//
	//     subnam, info, uplo, n, nb, imat
	//
	//
	//     subnam, info, uplo, n, nrhs, imat
	//
	//
	//     subnam, info, uplo, trans, diag, n, kd, nrhs, imat
	//
	//
	//     subnam, info, uplo, trans, diag, n, nrhs, imat
	//
	//
	//     subnam, info, uplo, trans, diag, NOrmin, n, imat
	//
	//
	//     subnam, info, uplo, trans, diag, NOrmin, n, kd, imat
	//
	//
	//     Unknown type
	//
	//
	//     What we do next
	//
	//
	return
	//
	//     End of Alaerh
	//
}
