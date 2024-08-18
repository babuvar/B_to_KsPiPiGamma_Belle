//               
//    cmlxerf.h   
//       ---- 
//    $Id: cmlxerf.h 9932 2006-11-12 14:26:53Z katayama $ 
//                
//    $Log$
//    Revision 1.1  2002/09/05 01:27:05  katayama
//    New package tatami from Nakadaira/Sumisawa san
//
//    Revision 1.2  2002/02/21 18:58:42  nakadair
//    add header
// 
//                

#ifndef __CMLXERF__H__
#define __CMLXERF__H__
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  double recexp_(double* x, double* y);
  double imcexp_(double* x, double* y);

  double rewerf_(double* x, double* y);
  double imwerf_(double* x, double* y);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __CMLXERF__H__ */
