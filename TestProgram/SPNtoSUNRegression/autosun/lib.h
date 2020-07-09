#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <tuple>
#include <map>
#include <vector>
using namespace std;
#include "./style.h"
#include "./list.h"
#include "./complex.h"
#include "./sparse.h"
#include "./polynomial.h"
#include "./matrix.h"
#include "./antisymmetric.h"
#include "./sun.h"

double* convertSPNToSUNAlgebra( const double* spnAlgebraVector){
    /**
     * given an spn algebra vector, retuns the corresponding 
     * sun algebra vector 
     */

    smatrix *TSUN,*TSPN;
    group::init(group::N,group::TYPESPN,TSPN);
    group::init(group::N,group::TYPESUN,TSUN);

    using group::N;
    double* res = new double[N*N-1];
    smatrix tmp;

    for(int isu = 0; isu<N*N-1;++isu){
        res[isu] = 0;
        for(int isp = 0; isp<N*(N+1)/2;++isp){
            complex tmpd;
            tmp.mult(TSPN[isu],TSUN[isp]);
            trace(tmpd,tmp);
            tmpd *= 2 ;  // because the normalization of SUN generators is 
                         // Tr(T^2) = 0.5
            if(tmpd.im != 0 ){
                cout << "tmpd.im != 0 " << endl;
                exit(1);
            }
            res[isu] += tmpd.re * spnAlgebraVector[isp];

        }
    }
    delete[] TSUN;
    delete[] TSPN;
    return res;
}

bool isnotzero(double x){
    const double threshold = 1e-14;
    return fabs(x)>threshold;
}



string write_spntosun_algebra(){
    /**
     * given an spn algebra vector, retuns the corresponding 
     * sun algebra vector 
     */

    smatrix *TSUN,*TSPN;
    group::init(group::N,group::TYPESPN,TSPN);
    group::init(group::N,group::TYPESUN,TSUN);


    stringstream ss;
    ss << "#define _spntosun_algebra(sunav,spnav) \\\n";

    using group::N;
    smatrix tmp;

    for(int isu = 0; isu<N*N-1;++isu){
        ss << "\tsunav[" << isu << "] =" ;
        for(int isp = 0; isp<N*(N+1)/2;++isp){
            complex tmpd;
            tmp.mult(TSUN[isu],TSPN[isp]);
            trace(tmpd,tmp);
            tmpd *= 2 ;  // because the normalization of SUN generators is 
                         // Tr(T^2) = 0.5
            if(isnotzero(tmpd.im)){
                cout << "tmpd.im != 0 " << endl;
                exit(1);
            }
            if(isnotzero(tmpd.re)){
                ss << " +("<< setprecision(15) << tmpd.re << "*spnav[" << isp << "])";
            }
        }
        ss << ";\\\n";
    }
    ss << "\n\n";
    delete[] TSUN;
    delete[] TSPN;
    return ss.str();
}
string write_suntospn_algebra(){
    /**
     * given an spn algebra vector, retuns the corresponding 
     * sun algebra vector 
     */

    smatrix *TSUN,*TSPN;
    group::init(group::N,group::TYPESPN,TSPN);
    group::init(group::N,group::TYPESUN,TSUN);


    stringstream ss;
    ss << "#define _suntospn_algebra(spnav,sunav) \\\n";

    using group::N;
    smatrix tmp;

    for(int isp = 0; isp<N*(N+1)/2;++isp){
        ss << "\tspnav[" << isp << "] =" ;
        for(int isu = 0; isu<N*N-1;++isu){
            complex tmpd;
            tmp.mult(TSPN[isp],TSUN[isu]);
            trace(tmpd,tmp);
            tmpd *= 2 ;  // because the normalization of SPN generators is 
                         // Tr(T^2) = 0.5
            if(isnotzero(tmpd.im)){
                cout << "tmpd.im != 0 " << endl;
                exit(1);
            }
            if(isnotzero(tmpd.re)){
                ss << " +("<< setprecision(15) << tmpd.re << "*sunav[" << isu << "])";
            }
        }
        ss << ";\\\n";
    }
    ss << "\n\n";
    delete[] TSUN;
    delete[] TSPN;
    return ss.str();
}

string write_conversions_base_tensors2(
        smatrix* TSRC, int source_group_algebra_dim, string src
        ,
        smatrix* TDST, int dest_group_algebra_dim, string dst 
        ,
        string reptype, 
        double basenorm
        ){
    using Entry =  tuple<int,double>;
    using Row = vector<Entry>;
    map<int,Row> dst_from_src;

    // computing scalar products between basis tensors  
    // in the spn and the sun adjoint representations
    // (i.e., the generators in the fundamental representation)
    for(int idst = 0; idst<dest_group_algebra_dim;++idst){
        Row row;
        for(int isrc = 0; isrc<source_group_algebra_dim;++isrc){
            complex tmpd;
            // scalar product between generators
            smatrix tmp;
            tmp.mult(TDST[idst],TSRC[isrc]);
            trace(tmpd,tmp);
            tmpd *= (double)1/basenorm ; 
            if(isnotzero(tmpd.im)){
                cout << "tmpd.im != 0 " << endl;
                exit(1);
            }
            if(isnotzero(tmpd.re)){
                row.push_back(Entry(isrc,tmpd.re));
            }
        }
        dst_from_src[idst] = row;
    }

    // now, using the transformation matrix computed for the 
    // fundamental representation, compute the transformation matrix
    // for a tensor having 2 indices in the adjoint of spn/sun
    // from its representation in the basis of the tensor product of spn/sun 
    // adjoint with itself to the basis of the tensor product of 
    // sun/spn adjoint with itself.
    
    auto dstIndex = [dest_group_algebra_dim](int i, int j){
        return i*dest_group_algebra_dim + j;
    };
    auto srcIndex = [source_group_algebra_dim](int i, int j){
        return i*source_group_algebra_dim + j;
    };

    stringstream ss;
    // For the adjoint:
    // src = spn, dst = sun:
    // "#define _spntnosun_adj(sunadjm,spnadjm) \\\n";
    // src = sun, dst = spn:
    // "#define _suntnospn_adj(spnadjm,sunadjm) \\\n";
    ss << "#define _" << src <<  "to" << dst << "_"<< reptype << "("
                      << dst << reptype << "m,"<< src << reptype << "m) \\\n";

    for(int isu1 = 0; isu1 < dest_group_algebra_dim; ++isu1)
        for(int isu2 = 0; isu2 < dest_group_algebra_dim; ++isu2){
            int rowIndex = dstIndex(isu1,isu2);
            ss << "\t"<< dst << reptype << "m[" << rowIndex << "]=";
            map<int,double> newRow; // 
            auto row1 = dst_from_src.at(isu1);
            auto row2 = dst_from_src.at(isu2);
            for(auto e1 : row1) for(auto e2 : row2){
                int c1,c2;
                double w1,w2;
                tie(c1,w1) = e1;
                tie(c2,w2) = e2;
                int colIndex = srcIndex(c1,c2);
                newRow[colIndex] += w1*w2; 
            }
            for(auto entry: newRow){
                int idx;
                double w;
                tie(idx,w) = entry;
                ss << "+(" << setprecision(15) << w << "*" << src << reptype << "m["
                    << idx << "])";
            }
        ss << ";\\\n";
        }
    ss << "\n\n";
    delete[] TSRC;
    delete[] TDST;
    return ss.str();
}

string write_spntosun_adj(){
    /**
     * Writes the macro that, given an spn algebra vector in the adjoint 
     * representation, expressed on the basis of the algebra of spn, return 
     * the components of it on the basis of the sun algebra vector in the 
     * adjoint representation.
     * The bases for the *fundamental representations* are created by 
     * group::init().
     */


    smatrix *TSUN,*TSPN;
    // create basis for SPN in the fundamental
    group::init(group::N,group::TYPESPN,TSPN);
    // create basis for SUN in the fundamental
    group::init(group::N,group::TYPESUN,TSUN);

    using group::N;
    int sunAlgebraCount = N*N-1;
    int spnAlgebraCount = N*(N+1)/2;

    return write_conversions_base_tensors2(
            TSPN,spnAlgebraCount,string("spn") // src
            ,
            TSUN,sunAlgebraCount,string("sun") // dst
            ,
            string("adj"), // reptype
            0.5 // the normalisation of all the generators
            );

}

string write_suntospn_adj(){
    /**
     * Writes the macro that, given an sun algebra vector in the adjoint 
     * representation, expressed on the basis of the algebra of sun, return 
     * the components of it on the basis of the spn algebra vector in the 
     * adjoint representation.
     * The bases for the *fundamental representations* are created by 
     * group::init().
     *
     * Notice that the space of algebra vectors for spn is a subspace 
     * of the space of algebra vectors for sun (so, the sun -> spn conversion
     * is actually a projection), and that both are a subspace of GL(N).
     */


    smatrix *TSUN,*TSPN;
    // create basis for SPN in the fundamental
    group::init(group::N,group::TYPESPN,TSPN);
    // create basis for SUN in the fundamental
    group::init(group::N,group::TYPESUN,TSUN);

    using group::N;
    int sunAlgebraCount = N*N-1;
    int spnAlgebraCount = N*(N+1)/2;

    return write_conversions_base_tensors2(
            TSUN,sunAlgebraCount,string("sun") // src
            ,                                
            TSPN,spnAlgebraCount,string("spn") // dst
            , 
            string("adj"), // reptype
            0.5 // the normalisation of all the generators
            );

}

string write_spntosun_asym(){
    /**
     * Writes the macro that, given a 2-index tensor in the antysimmetric 
     * representation of spn, returns the components of it on the basis of the 
     * antisymmetric representation of sun.
     * The bases for the *antisymmetric representations* are created by 
     * get_asymtensors_base*().
     */


    smatrix *ASYMSUN,*ASYMSPN;
    // create basis for SPN in the fundamental
    ASYMSPN = get_asymtensors_base_spn();
    // create basis for SUN in the fundamental
    ASYMSUN = get_asymtensors_base_sun();

    using group::N;
    int sunDim = N*(N-1)/2;
    int spnDim = sunDim -1;

    return write_conversions_base_tensors2(
            ASYMSPN,spnDim,string("spn") // src
            ,
            ASYMSUN,sunDim,string("sun") // dst
            , 
            string("asym"), // reptype
            1 // the normalisation of the base tensors
            );

}


string write_suntospn_asym(){
    /**
     * Writes the macro that, given a 2-index tensor in the antysimmetric 
     * representation of sun, returns the components of it on the basis of the 
     * antisymmetric representation of spn.
     * The bases for the *antisymmetric representations* are created by 
     * get_asymtensors_base*().
     *
     * Notice that the space of antisymmetric tensors for spn is a subspace 
     * of the space of antisymmetric tensors for sun (so, the sun -> spn 
     * conversion is actually a projection), and that both are a subspace of 
     * GL(N).
     */


    smatrix *ASYMSUN,*ASYMSPN;
    // create basis for SPN in the fundamental
    ASYMSPN = get_asymtensors_base_spn();
    // create basis for SUN in the fundamental
    ASYMSUN = get_asymtensors_base_sun();

    using group::N;
    int sunDim = N*(N-1)/2;
    int spnDim = sunDim -1;

    return write_conversions_base_tensors2(
            ASYMSUN,sunDim,string("sun") // src
            ,
            ASYMSPN,spnDim,string("spn") // dst
            , 
            string("asym"), // reptype
            1 // the normalisation of the base tensors
            );

}

void write_acf(string progname, string acfname, int N){ // algebra converter file
	group::N = N;
    ofstream acf(acfname.c_str()); 
    acf << "/** This file is produced automatically by " << progname << " */" << endl; 
    acf << "#define N " << N << endl << endl;
    acf << write_spntosun_adj();
    acf << write_suntospn_adj();
    acf << write_spntosun_asym();
    acf << write_suntospn_asym();
    acf << write_spntosun_algebra();
    acf << write_suntospn_algebra();
    acf.close();
}

void write_mtp_adj(string progname, string acfname, int N, string mtpname){
    ofstream mtp(mtpname.c_str()); // macro test program

    mtp << "/** This file is produced automatically by " << progname << " */" << endl; 
    mtp << "#include \"../colors.h\"" << endl;
    mtp << "#include <iostream>" << endl;
    mtp << "#include <iomanip>" << endl;
    mtp << "#include <cstdlib>" << endl;
    mtp << "#include <cmath>" << endl;
    mtp << "#include \"./" << acfname << "\"" << endl;
    mtp << "using namespace std;" << endl;

    mtp << "bool isnotzero(double x){" << endl;
    mtp << "    const double threshold = 1e-14;" << endl;
    mtp << "    return fabs(x)>threshold;" << endl;
    mtp << "}" << endl;
    mtp << "int main(){" << endl;
    int sunAlgebraCount = N*N-1;
    int spnAlgebraCount = N*(N+1)/2;

    mtp << "  {" << endl;
    mtp << "      bool ok = true;" << endl;
    mtp << "      double spnAlgebraVector["<<spnAlgebraCount<<"]; " << endl;
    mtp << "      double spnAlgebraVectorCheck["<<spnAlgebraCount<<"]; " << endl;
    mtp << "      double sunAlgebraVector["<<sunAlgebraCount<<"]; " << endl;
    mtp << "      for(int i = 0; i < "<<spnAlgebraCount<<"; ++i) " << endl;
    mtp << "          spnAlgebraVector[i] = (double)random()/RAND_MAX;" << endl;
    mtp << "      _spntosun_algebra(sunAlgebraVector,spnAlgebraVector);" << endl;
    mtp << "      _suntospn_algebra(spnAlgebraVectorCheck,sunAlgebraVector);" << endl;
    mtp << "      for(int i = 0; i < "<<spnAlgebraCount<<"; ++i){ " << endl;
    mtp << "          double check = spnAlgebraVectorCheck[i]-spnAlgebraVector[i];" << endl;
    mtp << "          if(isnotzero(check)){" << endl;
    mtp << "              cout << BOLDRED ;" << endl;
    mtp << "              cout << \"Problem With component \" << setw(4) << i << \" \" ;" << endl;
    mtp << "              cout << spnAlgebraVectorCheck[i] << \" vs  \" << spnAlgebraVector[i];" << endl;
    mtp << "              cout << endl;" << endl;
    mtp << "              cout << RESET;" << endl;
    mtp << "              ok=false;" << endl;
    mtp << "          };" << endl;
    mtp << "      };" << endl;
    mtp << "      if(ok) cout << \"All seems\" BOLDGREEN \" ok  \" RESET \"for the algebra.\" << endl; " << endl;
    mtp << "      else exit(1);" << endl;
    mtp << "  }" << endl;
    int spnAlgebraCountSq = spnAlgebraCount*spnAlgebraCount;
    int sunAlgebraCountSq = sunAlgebraCount*sunAlgebraCount;
    mtp << "  {" << endl;
    mtp << "      bool ok = true;" << endl;
    mtp << "      double spnAdjMatrix["<<spnAlgebraCountSq<<"]; " << endl;
    mtp << "      double spnAdjMatrixCheck["<<spnAlgebraCountSq <<"]; " << endl;
    mtp << "      double sunAdjMatrix["<<sunAlgebraCountSq <<"]; " << endl;
    mtp << "      for(int i = 0; i < "<<spnAlgebraCountSq<<"; ++i) " << endl;
    mtp << "          spnAdjMatrix[i] = (double)random()/RAND_MAX;" << endl;
    mtp << "      _spntosun_adj(sunAdjMatrix,spnAdjMatrix);" << endl;
    mtp << "      _suntospn_adj(spnAdjMatrixCheck,sunAdjMatrix);" << endl;
    mtp << "      for(int i = 0; i < "<<spnAlgebraCountSq<<"; ++i){ " << endl;
    mtp << "          double check = spnAdjMatrixCheck[i]-spnAdjMatrix[i];" << endl;
    mtp << "          if(isnotzero(check)){" << endl;
    mtp << "              cout << \"Problem With component \" << setw(4) << i << \" \" ;" << endl;
    mtp << "              cout << spnAdjMatrixCheck[i] << \" vs  \" << spnAdjMatrix[i];" << endl;
    mtp << "              cout << endl;" << endl;
    mtp << "              ok=false;" << endl;
    mtp << "          };" << endl;
    mtp << "      };" << endl;
    mtp << "      if(ok) cout << \"All seems\" BOLDGREEN \" ok  \" RESET \"for the adjmatrix.\" << endl; " << endl;
    mtp << "      else exit(1);" << endl;
    mtp << "  }" << endl;
    mtp << "  return 0; " << endl;
    mtp << "}" << endl;
    mtp.close();

}

void write_mtp_asym(string progname, string acfname, int N, string mtpname){
    ofstream mtp(mtpname.c_str()); // macro test program

    mtp << "/** This file is produced automatically by " << progname << " */" << endl; 
    mtp << "#include \"../colors.h\"" << endl;
    mtp << "#include <iostream>" << endl;
    mtp << "#include <iomanip>" << endl;
    mtp << "#include <cstdlib>" << endl;
    mtp << "#include <cmath>" << endl;
    mtp << "#include \"./" << acfname << "\"" << endl;
    mtp << "using namespace std;" << endl;

    mtp << "bool isnotzero(double x){" << endl;
    mtp << "    const double threshold = 1e-14;" << endl;
    mtp << "    return fabs(x)>threshold;" << endl;
    mtp << "}" << endl;
    mtp << "int main(){" << endl;
    int sunReprDim = N*(N-1)/2;
    int spnReprDim = sunReprDim - 1;
    int spnReprDimSq = spnReprDim*spnReprDim;
    int sunReprDimSq = sunReprDim*sunReprDim;
    mtp << "  {" << endl;
    mtp << "      bool ok = true;" << endl;
    mtp << "      double spnAsymMatrix["<<spnReprDimSq<<"]; " << endl;
    mtp << "      double spnAsymMatrixCheck["<<spnReprDimSq <<"]; " << endl;
    mtp << "      double sunAsymMatrix["<<sunReprDimSq <<"]; " << endl;
    mtp << "      for(int i = 0; i < "<<spnReprDimSq<<"; ++i) " << endl;
    mtp << "          spnAsymMatrix[i] = (double)random()/RAND_MAX;" << endl;
    mtp << "      _spntosun_asym(sunAsymMatrix,spnAsymMatrix);" << endl;
    mtp << "      _suntospn_asym(spnAsymMatrixCheck,sunAsymMatrix);" << endl;
    mtp << "      for(int i = 0; i < "<<spnReprDimSq<<"; ++i){ " << endl;
    mtp << "          double check = spnAsymMatrixCheck[i]-spnAsymMatrix[i];" << endl;
    mtp << "          if(isnotzero(check)){" << endl;
    mtp << "              cout << BOLDRED ;" << endl;
    mtp << "              cout << \"Problem With component \" << setw(4) << i << \" \" ;" << endl;
    mtp << "              cout << spnAsymMatrixCheck[i] << \" vs  \" << spnAsymMatrix[i];" << endl;
    mtp << "              cout << RESET;" << endl;
    mtp << "              cout << endl;" << endl;
    mtp << "              ok=false;" << endl;
    mtp << "          };" << endl;
    mtp << "      };" << endl;
    mtp << "      if(ok) cout << \"All seems\" BOLDGREEN \" ok  \" RESET \"for the asym matrix.\" << endl; " << endl;
    mtp << "      else exit(1);" << endl;
    mtp << "  }" << endl;
    mtp << "  return 0; " << endl;
    mtp << "}" << endl;
    mtp.close();

}

