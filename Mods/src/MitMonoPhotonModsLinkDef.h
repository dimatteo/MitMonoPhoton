// $Id: MitMonoPhotonModsLinkDef.h,v 1.2 2013/06/23 14:31:17 ceballos Exp $

#ifndef MITMONOPHOTON_MODS_LINKDEF_H
#define MITMONOPHOTON_MODS_LINKDEF_H
#include "MitMonoPhoton/Mods/interface/MonoPhotonTreeWriter.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MonoPhotonTreeWriter+;
#endif
