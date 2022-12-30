#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>

         template<typename T, std::size_t N>
         std::size_t length_of( T const (&)[N] )
         {
             return N;
         }
 
         void initializeCartesianExponentsMoldenFormat(std::vector<std::vector<int> > & sCartExp,
                                                       std::vector<std::vector<int> > & pCartExp,
                                                       std::vector<std::vector<int> > & dCartExp,
                                                       std::vector<std::vector<int> > & fCartExp,
                                                       std::vector<std::vector<int> > & gCartExp)
         {
 
           const char * pCartExpChars[] = {"x", "y", "z"};
           const char * dCartExpChars[] = {"xx", "yy", "zz", "xy", "xz", "yz"};
           const char * fCartExpChars[] = {"xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"};
           const char * gCartExpChars[] = {"xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy",
                                             "xxyy", "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"};
 
           int pSize = length_of(pCartExpChars);
           int dSize = length_of(dCartExpChars);
           int fSize = length_of(fCartExpChars);
           int gSize = length_of(gCartExpChars);
 
           assert(pSize == 3);
           assert(dSize == 6);
           assert(fSize == 10);
           assert(gSize == 15);
 
           std::vector<std::string> pCartExpStrings(pCartExpChars, pCartExpChars + pSize);
           std::vector<std::string> dCartExpStrings(dCartExpChars, dCartExpChars + dSize);
           std::vector<std::string> fCartExpStrings(fCartExpChars, fCartExpChars + fSize);
           std::vector<std::string> gCartExpStrings(gCartExpChars, gCartExpChars + gSize);
 
           sCartExp.resize(1, std::vector<int>(3,0));
           pCartExp.resize(pSize, std::vector<int>(3));
           dCartExp.resize(dSize, std::vector<int>(3));
           fCartExp.resize(fSize, std::vector<int>(3));
           gCartExp.resize(gSize, std::vector<int>(3));
 
           for(unsigned int i = 0; i < pSize; ++i)
           {
             pCartExp[i][0] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'x');
             pCartExp[i][1] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'y');
             pCartExp[i][2] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'z');

           }
              for(unsigned int i = 0; i < dSize; ++i)
           {
             dCartExp[i][0] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'x');
             dCartExp[i][1] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'y');
             dCartExp[i][2] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'z');
           }
 
           for(unsigned int i = 0; i < fSize; ++i)
           {
             fCartExp[i][0] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'x');
             fCartExp[i][1] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'y');
             fCartExp[i][2] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'z');
           }
 
           for(unsigned int i = 0; i < gSize; ++i)
           {
             gCartExp[i][0] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'x');
             gCartExp[i][1] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'y');
             gCartExp[i][2] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'z');
           }
 
           //delete pCartExpChars;
           //delete dCartExpChars;
           //delete fCartExpChars;
           //delete gCartExpChars;

        }

        int main()
        {

            std::vector<std::vector<int> >  sCartExp;
            std::vector<std::vector<int> >  pCartExp;
            std::vector<std::vector<int> >  dCartExp;
            std::vector<std::vector<int> >  fCartExp;
            std::vector<std::vector<int> >  gCartExp;
            initializeCartesianExponentsMoldenFormat(sCartExp,
                                                     pCartExp,
                                                     dCartExp,
                                                     fCartExp,
                                                     gCartExp);   
        }
