// $Id: MitMonoPhotonUtilsLinkDef.h,v 1.1 2013/08/22 22:09:43 dimatteo Exp $

#ifndef MITMONOPHOTON_UTILS_LINKDEF_H
#define MITMONOPHOTON_UTILS_LINKDEF_H
#include "MitMonoPhoton/Utils/interface/TreeReducer.h"
#include "MitMonoPhoton/Utils/interface/TreeReducerFake.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::TreeReducer+;
#pragma link C++ class mithep::TreeReducerFake+;
#endif
