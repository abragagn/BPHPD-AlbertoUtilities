#ifndef TopDecayMode_h
#define TopDecayMode_h

class TopDecayMode {

 public:

  enum type { twl  = 11, twtl,
              tbl  = 21, tbtl, tbbl, tbcl = 25, tbol, tbxl,
              twcl = 31, twxl,
              bl   = 41, btl , bbl , bcl  = 45 , bol,  bxl,
              cl   = 51, xl  ,
              other = 99,
              tb    = 111, twq,
              xb    = 121,
              hadronic = 199,
              unknown  = 999 };

};

#endif

