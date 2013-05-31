// $Id: MitPhysicsModsLinkDef.h,v 1.38 2012/11/14 20:05:37 sixie Exp $

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
#pragma link C++ class mithep::MonoPhotonEvent+;
#endif
