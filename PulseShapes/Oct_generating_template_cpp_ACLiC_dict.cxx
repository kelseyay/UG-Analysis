//
// File generated by /usr/bin/rootcint at Wed May 30 22:51:28 2018

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dIafsdIcerndOchdIuserdIkdIkyeedIUGmIAnalysisdIPulseShapesdIOct_generating_template_cpp_ACLiC_dict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "Oct_generating_template_cpp_ACLiC_dict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

/********************************************************
* /afs/cern.ch/user/k/kyee/UG-Analysis/PulseShapes/Oct_generating_template_cpp_ACLiC_dict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableOct_generating_template_cpp_ACLiC_dict();

extern "C" void G__set_cpp_environmentOct_generating_template_cpp_ACLiC_dict() {
  G__cpp_reset_tagtableOct_generating_template_cpp_ACLiC_dict();
}
#include <new>
extern "C" int G__cpp_dllrevOct_generating_template_cpp_ACLiC_dict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* Setting up global function */
static int G__Oct_generating_template_cpp_ACLiC_dict__0_2858(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      Oct_generating_template();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/*********************************************************
* Member function Stub
*********************************************************/

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncOct_generating_template_cpp_ACLiC_dict {
 public:
  G__Sizep2memfuncOct_generating_template_cpp_ACLiC_dict(): p(&G__Sizep2memfuncOct_generating_template_cpp_ACLiC_dict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncOct_generating_template_cpp_ACLiC_dict::*p)();
};

size_t G__get_sizep2memfuncOct_generating_template_cpp_ACLiC_dict()
{
  G__Sizep2memfuncOct_generating_template_cpp_ACLiC_dict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceOct_generating_template_cpp_ACLiC_dict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableOct_generating_template_cpp_ACLiC_dict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Float_t>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Double_t>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("pair<UInt_t,Int_t>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_pairlEunsignedsPintcOintgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<std::pair<UInt_t,Int_t> >",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<pair<UInt_t,Int_t> >",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<pair<unsigned int,int> >",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<std::string,TObjArray*>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*,less<string> >",117,G__get_linked_tagnum(&G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */
extern "C" void G__cpp_setup_memvarOct_generating_template_cpp_ACLiC_dict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/

/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncOct_generating_template_cpp_ACLiC_dict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {
}

static void G__cpp_setup_global4() {
}

static void G__cpp_setup_global5() {
   G__memvar_setup((void*)(&nsamples),105,0,1,-1,-1,-1,1,"nsamples=",0,(char*)NULL);
   G__memvar_setup((void*)(&nboards),105,0,1,-1,-1,-1,1,"nboards=",0,(char*)NULL);
   G__memvar_setup((void*)(&nchannels_per_board),105,0,1,-1,-1,-1,1,"nchannels_per_board=",0,(char*)NULL);
   G__memvar_setup((void*)(&nchannels),105,0,1,-1,-1,-1,1,"nchannels=",0,(char*)NULL);
   G__memvar_setup((void*)(&impinged_channel),105,0,1,-1,-1,-1,1,"impinged_channel=",0,(char*)NULL);
   G__memvar_setup((void*)(&beam_energy),102,0,1,-1,-1,-1,1,"beam_energy=",0,(char*)NULL);
   G__memvar_setup((void*)(&amp_max),102,0,0,-1,-1,-1,1,"amp_max[5]=",0,(char*)NULL);
   G__memvar_setup((void*)(&time_max),102,0,0,-1,-1,-1,1,"time_max[5]=",0,(char*)NULL);
   G__memvar_setup((void*)(&WF_time),102,0,0,-1,-1,-1,1,"WF_time[750]=",0,(char*)NULL);
   G__memvar_setup((void*)(&WF_val),102,0,0,-1,-1,-1,1,"WF_val[750]=",0,(char*)NULL);
   G__memvar_setup((void*)(&X),102,0,0,-1,-1,-1,1,"X[2]=",0,(char*)NULL);
   G__memvar_setup((void*)(&Y),102,0,0,-1,-1,-1,1,"Y[2]=",0,(char*)NULL);
   G__memvar_setup((void*)(&hodoX),102,0,0,-1,-1,-1,1,"hodoX[2]=",0,(char*)NULL);
   G__memvar_setup((void*)(&hodoY),102,0,0,-1,-1,-1,1,"hodoY[2]=",0,(char*)NULL);

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalOct_generating_template_cpp_ACLiC_dict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
  G__cpp_setup_global5();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {
}

static void G__cpp_setup_func23() {
}

static void G__cpp_setup_func24() {
}

static void G__cpp_setup_func25() {
}

static void G__cpp_setup_func26() {
}

static void G__cpp_setup_func27() {
}

static void G__cpp_setup_func28() {
   G__memfunc_setup("Oct_generating_template", 2404, G__Oct_generating_template_cpp_ACLiC_dict__0_2858, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL
, (void*) NULL, 0);

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcOct_generating_template_cpp_ACLiC_dict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
  G__cpp_setup_func23();
  G__cpp_setup_func24();
  G__cpp_setup_func25();
  G__cpp_setup_func26();
  G__cpp_setup_func27();
  G__cpp_setup_func28();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEdoublegR = { "TMatrixTBase<double>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEfloatgR = { "TMatrixTBase<float>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_pairlEunsignedsPintcOintgR = { "pair<unsigned int,int>" , 115 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR = { "vector<pair<unsigned int,int>,allocator<pair<unsigned int,int> > >" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<pair<unsigned int,int>,allocator<pair<unsigned int,int> > >::iterator>" , 99 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR = { "map<string,TObjArray*,less<string>,allocator<pair<const string,TObjArray*> > >" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableOct_generating_template_cpp_ACLiC_dict() {
  G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEdoublegR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEfloatgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_pairlEunsignedsPintcOintgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableOct_generating_template_cpp_ACLiC_dict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEdoublegR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_TMatrixTBaselEfloatgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEdoublegR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_pairlEunsignedsPintcOintgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_vectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_reverse_iteratorlEvectorlEpairlEunsignedsPintcOintgRcOallocatorlEpairlEunsignedsPintcOintgRsPgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__Oct_generating_template_cpp_ACLiC_dictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR);
}
extern "C" void G__cpp_setupOct_generating_template_cpp_ACLiC_dict(void) {
  G__check_setup_version(30051515,"G__cpp_setupOct_generating_template_cpp_ACLiC_dict()");
  G__set_cpp_environmentOct_generating_template_cpp_ACLiC_dict();
  G__cpp_setup_tagtableOct_generating_template_cpp_ACLiC_dict();

  G__cpp_setup_inheritanceOct_generating_template_cpp_ACLiC_dict();

  G__cpp_setup_typetableOct_generating_template_cpp_ACLiC_dict();

  G__cpp_setup_memvarOct_generating_template_cpp_ACLiC_dict();

  G__cpp_setup_memfuncOct_generating_template_cpp_ACLiC_dict();
  G__cpp_setup_globalOct_generating_template_cpp_ACLiC_dict();
  G__cpp_setup_funcOct_generating_template_cpp_ACLiC_dict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncOct_generating_template_cpp_ACLiC_dict();
  return;
}
class G__cpp_setup_initOct_generating_template_cpp_ACLiC_dict {
  public:
    G__cpp_setup_initOct_generating_template_cpp_ACLiC_dict() { G__add_setup_func("Oct_generating_template_cpp_ACLiC_dict",(G__incsetup)(&G__cpp_setupOct_generating_template_cpp_ACLiC_dict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initOct_generating_template_cpp_ACLiC_dict() { G__remove_setup_func("Oct_generating_template_cpp_ACLiC_dict"); }
};
G__cpp_setup_initOct_generating_template_cpp_ACLiC_dict G__cpp_setup_initializerOct_generating_template_cpp_ACLiC_dict;
