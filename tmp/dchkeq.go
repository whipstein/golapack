package goblas

import (
	"math"
)

// Dchkeq tests Dgeequ, Dgbequ, Dpoequ, Dppequ and Dpbequ
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkeq( thresh, nout)
//
//       .. Scalar Arguments ..
//       intEGER            nout
//       DOUBLE PRECISION   thresh
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkeq tests Dgeequ, Dgbequ, Dpoequ, Dppequ and Dpbequ
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] thresh
// \verbatim
//          thresh is DOUBLE PRECISION
//          Threshold for testing routines. Should be between 2 and 10.
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number for output.
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
func Dchkeq(thresh *float64, nout *int) {
	zero := new(float64)
	one := new(float64)
	ten := new(float64)
	nsz := new(int)
	ns2b := new(int)
	nszp := new(int)
	npow := new(int)
	ok := new(bool)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	info := new(int)
	j := new(int)
	kl := new(int)
	ku := new(int)
	m := new(int)
	n := new(int)
	ccond := new(float64)
	eps := new(float64)
	norm := new(float64)
	ratio := new(float64)
	rcmax := new(float64)
	rcmin := new(float64)
	rcond := new(float64)
	a := func() *[][]float64 {
		arr := make([][]float64, 5)
		for u := 0; u < 5; u++ {
			arr[u] = make([]float64, 5)
		}
		return &arr
	}()
	ab := func() *[][]float64 {
		arr := make([][]float64, 0)
		for u := 0; u < 0; u++ {
			arr[u] = make([]float64, 5)
		}
		return &arr
	}()
	ap := func() *[]float64 {
		arr := make([]float64, 0)
		return &arr
	}()
	c := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	pow := func() *[]float64 {
		arr := make([]float64, 0)
		return &arr
	}()
	r := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	reslts := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	rpow := func() *[]float64 {
		arr := make([]float64, 0)
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
	(*zero) = 0.0
	(*one) = 1.0e+0
	(*ten) = 1.0e1
	(*nsz) = 5
	(*ns2b) = 3*(*nsz) - 2
	(*nszp) = ((*nsz) * ((*nsz) + 1)) / 2
	(*npow) = 2*(*nsz) + 1
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	(*path)[0] = *func() *[]byte {y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y := []byte("EQ"); return &y }()
	//
	(*eps) = (*Dlamch(func() *byte {y := byte('P'); return &y }()))
	for (*i) = 1; (*i) <= 5; (*i)++ {
		(*reslts)[(*i)-1] = (*zero)
		//Label10:
	}
	for (*i) = 1; (*i) <= (*npow); (*i)++ {
		(*pow)[(*i)-1] = math.Pow((*ten), ((*i) - 1))
		(*rpow)[(*i)-1] = (*one) / (*pow)[(*i)-1]
		//Label20:
	}
	//
	//     Test Dgeequ
	//
	for (*n) = 0; (*n) <= (*nsz); (*n)++ {
		for (*m) = 0; (*m) <= (*nsz); (*m)++ {
			//
			for (*j) = 1; (*j) <= (*nsz); (*j)++ {
				for (*i) = 1; (*i) <= (*nsz); (*i)++ {
					if (*i) <= (*m) && (*j) <= (*n) {
						(*a)[(*i)-1][(*j)-1] = (*pow)[(*i)+(*j)+0] * math.Pow((-1), ((*i)+(*j)))
					} else {
						(*a)[(*i)-1][(*j)-1] = (*zero)
					}
					//Label30:
				}
				//Label40:
			}
			//
			Dgeequ(m, n, a, nsz, r, c, rcond, ccond, norm, info)
			//
			if (*info) != 0 {
				(*reslts)[0] = (*one)
			} else {
				if (*n) != 0 && (*m) != 0 {
					(*reslts)[0] = (MAX(((*reslts)[0]), ABS(((*rcond)-(*rpow)[(*m)-1])/(*rpow)[(*m)-1])))
					(*reslts)[0] = (MAX(((*reslts)[0]), ABS(((*ccond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
					(*reslts)[0] = (MAX(((*reslts)[0]), ABS(((*norm)-(*pow)[(*n)+(*m)+0])/(*pow)[(*n)+(*m)+0])))
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*reslts)[0] = (MAX(((*reslts)[0]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+(*n)+0])/(*rpow)[(*i)+(*n)+0])))
						//Label50:
					}
					for (*j) = 1; (*j) <= (*n); (*j)++ {
						(*reslts)[0] = (MAX(((*reslts)[0]), ABS(((*c)[(*j)-1]-(*pow)[(*n)-(*j)+0])/(*pow)[(*n)-(*j)+0])))
						//Label60:
					}
				}
			}
			//
			//Label70:
		}
		//Label80:
	}
	//
	//     Test with zero rows and columns
	//
	for (*j) = 1; (*j) <= (*nsz); (*j)++ {
		(*a)[MAX((*nsz)-1, 1)-1][(*j)-1] = (*zero)
		//Label90:
	}
	Dgeequ(nsz, nsz, a, nsz, r, c, rcond, ccond, norm, info)
	if (*info) != (MAX((*nsz)-1, 1)) {
		(*reslts)[0] = (*one)
	}
	//
	for (*j) = 1; (*j) <= (*nsz); (*j)++ {
		(*a)[MAX((*nsz)-1, 1)-1][(*j)-1] = (*one)
		//Label100:
	}
	for (*i) = 1; (*i) <= (*nsz); (*i)++ {
		(*a)[(*i)-1][MAX((*nsz)-1, 1)-1] = (*zero)
		//Label110:
	}
	Dgeequ(nsz, nsz, a, nsz, r, c, rcond, ccond, norm, info)
	if (*info) != (*nsz)+MAX((*nsz)-1, 1) {
		(*reslts)[0] = (*one)
	}
	(*reslts)[0] = (*reslts)[0] / (*eps)
	//
	//     Test Dgbequ
	//
	for (*n) = 0; (*n) <= (*nsz); (*n)++ {
		for (*m) = 0; (*m) <= (*nsz); (*m)++ {
			for (*kl) = 0; (*kl) <= (MAX((*m)-1, 0)); (*kl)++ {
				for (*ku) = 0; (*ku) <= (MAX((*n)-1, 0)); (*ku)++ {
					//
					for (*j) = 1; (*j) <= (*nsz); (*j)++ {
						for (*i) = 1; (*i) <= (*ns2b); (*i)++ {
							(*ab)[(*i)-1][(*j)-1] = (*zero)
							//Label120:
						}
						//Label130:
					}
					for (*j) = 1; (*j) <= (*n); (*j)++ {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							if (*i) <= Min((*m), (*j)+(*kl)) && (*i) >= MAX(1, (*j)-(*ku)) && (*j) <= (*n) {
								(*ab)[(*ku)+1+(*i)-(*j)-1][(*j)-1] = (*pow)[(*i)+(*j)+0] * math.Pow((-1), ((*i)+(*j)))
							}
							//Label140:
						}
						//Label150:
					}
					//
					Dgbequ(m, n, kl, ku, ab, ns2b, r, c, rcond, ccond, norm, info)
					//
					if (*info) != 0 {
						if !(((*n)+(*kl) < (*m) && (*info) == (*n)+(*kl)+1) || ((*m)+(*ku) < (*n) && (*info) == 2*(*m)+(*ku)+1)) {
							(*reslts)[1] = (*one)
						}
					} else {
						if (*n) != 0 && (*m) != 0 {
							//
							(*rcmin) = (*r)[0]
							(*rcmax) = (*r)[0]
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*rcmin) = (Min((*rcmin), ((*r)[(*i)-1])))
								(*rcmax) = (MAX((*rcmax), ((*r)[(*i)-1])))
								//Label160:
							}
							(*ratio) = (*rcmin) / (*rcmax)
							(*reslts)[1] = (MAX(((*reslts)[1]), ABS(((*rcond)-(*ratio))/(*ratio))))
							//
							(*rcmin) = (*c)[0]
							(*rcmax) = (*c)[0]
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								(*rcmin) = (Min((*rcmin), ((*c)[(*j)-1])))
								(*rcmax) = (MAX((*rcmax), ((*c)[(*j)-1])))
								//Label170:
							}
							(*ratio) = (*rcmin) / (*rcmax)
							(*reslts)[1] = (MAX(((*reslts)[1]), ABS(((*ccond)-(*ratio))/(*ratio))))
							//
							(*reslts)[1] = (MAX(((*reslts)[1]), ABS(((*norm)-(*pow)[(*n)+(*m)+0])/(*pow)[(*n)+(*m)+0])))
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*rcmax) = (*zero)
								for (*j) = 1; (*j) <= (*n); (*j)++ {
									if (*i) <= (*j)+(*kl) && (*i) >= (*j)-(*ku) {
										(*ratio) = (ABS((*r)[(*i)-1] * (*pow)[(*i)+(*j)+0] * (*c)[(*j)-1]))
										(*rcmax) = (MAX((*rcmax), (*ratio)))
									}
									//Label180:
								}
								(*reslts)[1] = (MAX(((*reslts)[1]), ABS((*one)-(*rcmax))))
								//Label190:
							}
							//
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								(*rcmax) = (*zero)
								for (*i) = 1; (*i) <= (*m); (*i)++ {
									if (*i) <= (*j)+(*kl) && (*i) >= (*j)-(*ku) {
										(*ratio) = (ABS((*r)[(*i)-1] * (*pow)[(*i)+(*j)+0] * (*c)[(*j)-1]))
										(*rcmax) = (MAX((*rcmax), (*ratio)))
									}
									//Label200:
								}
								(*reslts)[1] = (MAX(((*reslts)[1]), ABS((*one)-(*rcmax))))
								//Label210:
							}
						}
					}
					//
					//Label220:
				}
				//Label230:
			}
			//Label240:
		}
		//Label250:
	}
	(*reslts)[1] = (*reslts)[1] / (*eps)
	//
	//     Test Dpoequ
	//
	for (*n) = 0; (*n) <= (*nsz); (*n)++ {
		//
		for (*i) = 1; (*i) <= (*nsz); (*i)++ {
			for (*j) = 1; (*j) <= (*nsz); (*j)++ {
				if (*i) <= (*n) && (*j) == (*i) {
					(*a)[(*i)-1][(*j)-1] = (*pow)[(*i)+(*j)+0] * math.Pow((-1), ((*i)+(*j)))
				} else {
					(*a)[(*i)-1][(*j)-1] = (*zero)
				}
				//Label260:
			}
			//Label270:
		}
		//
		Dpoequ(n, a, nsz, r, rcond, norm, info)
		//
		if (*info) != 0 {
			(*reslts)[2] = (*one)
		} else {
			if (*n) != 0 {
				(*reslts)[2] = (MAX(((*reslts)[2]), ABS(((*rcond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
				(*reslts)[2] = (MAX(((*reslts)[2]), ABS(((*norm)-(*pow)[2*(*n)+0])/(*pow)[2*(*n)+0])))
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*reslts)[2] = (MAX(((*reslts)[2]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+0])/(*rpow)[(*i)+0])))
					//Label280:
				}
			}
		}
		//Label290:
	}
	(*a)[MAX((*nsz)-1, 1)-1][MAX((*nsz)-1, 1)-1] = -(*one)
	Dpoequ(nsz, a, nsz, r, rcond, norm, info)
	if (*info) != (MAX((*nsz)-1, 1)) {
		(*reslts)[2] = (*one)
	}
	(*reslts)[2] = (*reslts)[2] / (*eps)
	//
	//     Test Dppequ
	//
	for (*n) = 0; (*n) <= (*nsz); (*n)++ {
		//
		//        Upper triangular packed storage
		//
		for (*i) = 1; (*i) <= ((*n)*((*n)+1))/2; (*i)++ {
			(*ap)[(*i)-1] = (*zero)
			//Label300:
		}
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ap)[((*i)*((*i)+1))/1] = (*pow)[2*(*i)+0]
			//Label310:
		}
		//
		Dppequ(func() *byte {y := byte('U'); return &y }(), n, ap, r, rcond, norm, info)
		//
		if (*info) != 0 {
			(*reslts)[3] = (*one)
		} else {
			if (*n) != 0 {
				(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*rcond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
				(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*norm)-(*pow)[2*(*n)+0])/(*pow)[2*(*n)+0])))
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+0])/(*rpow)[(*i)+0])))
					//Label320:
				}
			}
		}
		//
		//        Lower triangular packed storage
		//
		for (*i) = 1; (*i) <= ((*n)*((*n)+1))/2; (*i)++ {
			(*ap)[(*i)-1] = (*zero)
			//Label330:
		}
		(*j) = 1
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ap)[(*j)-1] = (*pow)[2*(*i)+0]
			(*j) = (*j) + ((*n) - (*i) + 1)
			//Label340:
		}
		//
		Dppequ(func() *byte {y := byte('L'); return &y }(), n, ap, r, rcond, norm, info)
		//
		if (*info) != 0 {
			(*reslts)[3] = (*one)
		} else {
			if (*n) != 0 {
				(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*rcond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
				(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*norm)-(*pow)[2*(*n)+0])/(*pow)[2*(*n)+0])))
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*reslts)[3] = (MAX(((*reslts)[3]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+0])/(*rpow)[(*i)+0])))
					//Label350:
				}
			}
		}
		//
		//Label360:
	}
	(*i) = ((*nsz)*((*nsz)+1))/2 - 2
	(*ap)[(*i)-1] = -(*one)
	Dppequ(func() *byte {y := byte('L'); return &y }(), nsz, ap, r, rcond, norm, info)
	if (*info) != (MAX((*nsz)-1, 1)) {
		(*reslts)[3] = (*one)
	}
	(*reslts)[3] = (*reslts)[3] / (*eps)
	//
	//     Test Dpbequ
	//
	for (*n) = 0; (*n) <= (*nsz); (*n)++ {
		for (*kl) = 0; (*kl) <= (MAX((*n)-1, 0)); (*kl)++ {
			//
			//           Test upper triangular storage
			//
			for (*j) = 1; (*j) <= (*nsz); (*j)++ {
				for (*i) = 1; (*i) <= (*ns2b); (*i)++ {
					(*ab)[(*i)-1][(*j)-1] = (*zero)
					//Label370:
				}
				//Label380:
			}
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*ab)[(*kl)+0][(*j)-1] = (*pow)[2*(*j)+0]
				//Label390:
			}
			//
			Dpbequ(func() *byte {y := byte('U'); return &y }(), n, kl, ab, ns2b, r, rcond, norm, info)
			//
			if (*info) != 0 {
				(*reslts)[4] = (*one)
			} else {
				if (*n) != 0 {
					(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*rcond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
					(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*norm)-(*pow)[2*(*n)+0])/(*pow)[2*(*n)+0])))
					for (*i) = 1; (*i) <= (*n); (*i)++ {
						(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+0])/(*rpow)[(*i)+0])))
						//Label400:
					}
				}
			}
			if (*n) != 0 {
				(*ab)[(*kl)+0][MAX((*n)-1, 1)-1] = -(*one)
				Dpbequ(func() *byte {y := byte('U'); return &y }(), n, kl, ab, ns2b, r, rcond, norm, info)
				if (*info) != (MAX((*n)-1, 1)) {
					(*reslts)[4] = (*one)
				}
			}
			//
			//           Test lower triangular storage
			//
			for (*j) = 1; (*j) <= (*nsz); (*j)++ {
				for (*i) = 1; (*i) <= (*ns2b); (*i)++ {
					(*ab)[(*i)-1][(*j)-1] = (*zero)
					//Label410:
				}
				//Label420:
			}
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*ab)[0][(*j)-1] = (*pow)[2*(*j)+0]
				//Label430:
			}
			//
			Dpbequ(func() *byte {y := byte('L'); return &y }(), n, kl, ab, ns2b, r, rcond, norm, info)
			//
			if (*info) != 0 {
				(*reslts)[4] = (*one)
			} else {
				if (*n) != 0 {
					(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*rcond)-(*rpow)[(*n)-1])/(*rpow)[(*n)-1])))
					(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*norm)-(*pow)[2*(*n)+0])/(*pow)[2*(*n)+0])))
					for (*i) = 1; (*i) <= (*n); (*i)++ {
						(*reslts)[4] = (MAX(((*reslts)[4]), ABS(((*r)[(*i)-1]-(*rpow)[(*i)+0])/(*rpow)[(*i)+0])))
						//Label440:
					}
				}
			}
			if (*n) != 0 {
				(*ab)[0][MAX((*n)-1, 1)-1] = -(*one)
				Dpbequ(func() *byte {y := byte('L'); return &y }(), n, kl, ab, ns2b, r, rcond, norm, info)
				if (*info) != (MAX((*n)-1, 1)) {
					(*reslts)[4] = (*one)
				}
			}
			//Label450:
		}
		//Label460:
	}
	(*reslts)[4] = (*reslts)[4] / (*eps)
	(*ok) = ((*reslts)[0] <= (*(thresh))) && ((*reslts)[1] <= (*(thresh))) && ((*reslts)[2] <= (*(thresh))) && ((*reslts)[3] <= (*(thresh))) && ((*reslts)[4] <= (*(thresh)))
	WRITE((*(nout)), *func() *[]byte {y := []byte(" %v\n"); return &y }())
	if *ok {
		WRITE((*(nout)), *func() *[]byte {y := []byte(" All tests for %3s routines passed the threshold\n"); return &y }(), (*path))
	} else {
		if (*reslts)[0] > (*(thresh)) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" Dgeequ failed test with value %10.3E exceeding threshold %10.3E\n")
				return &y
			}(), (*reslts)[0], (*(thresh)))
		}
		if (*reslts)[1] > (*(thresh)) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" Dgbequ failed test with value %10.3E exceeding threshold %10.3E\n")
				return &y
			}(), (*reslts)[1], (*(thresh)))
		}
		if (*reslts)[2] > (*(thresh)) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" Dpoequ failed test with value %10.3E exceeding threshold %10.3E\n")
				return &y
			}(), (*reslts)[2], (*(thresh)))
		}
		if (*reslts)[3] > (*(thresh)) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" Dppequ failed test with value %10.3E exceeding threshold %10.3E\n")
				return &y
			}(), (*reslts)[3], (*(thresh)))
		}
		if (*reslts)[4] > (*(thresh)) {
			WRITE((*(nout)), *func() *[]byte {
				y := []byte(" Dpbequ failed test with value %10.3E exceeding threshold %10.3E\n")
				return &y
			}(), (*reslts)[4], (*(thresh)))
		}
	}
	return
	//
	//     End of Dchkeq
	//
}
