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

enum CONV_DIRECTION { PTOU // sPn to sUn
                    , UTOP // sUn to sPn
                    };

string write_conversions_adj_tensors(CONV_DIRECTION dir){
    /**
     * Writes the macro that, given an spn/sun algebra vector in the adjoint 
     * representation, expressed on the basis of the algebra of spn/sun, return 
     * the components of it on the basis of the sun/spn algebra vector in the 
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
    smatrix tmp;

    using Entry =  tuple<int,double>;
    using Row = vector<Entry>;
    map<int,Row> dst_from_src;

    int dest_group_algebra_dim;
    int source_group_algebra_dim;
    string src,dst;
    {
        int sunAlgebraCount = N*N-1;
        int spnAlgebraCount = N*(N+1)/2;
        switch(dir){
            case PTOU:
                dest_group_algebra_dim = sunAlgebraCount;
                source_group_algebra_dim = spnAlgebraCount;
                dst = "sun";
                src = "spn";
                break;
            case UTOP:
                dest_group_algebra_dim = spnAlgebraCount;
                source_group_algebra_dim = sunAlgebraCount;
                dst = "spn";
                src = "sun";
                break;
            default:
                exit(1);
        }
    }
    // computing scalar products between generators
    // in the spn and the sun algebra - in the fundamental representation
    for(int idst = 0; idst<dest_group_algebra_dim;++idst){
        Row row;
        for(int isrc = 0; isrc<source_group_algebra_dim;++isrc){
            complex tmpd;
            // scalar product between generators
            switch(dir){
                case PTOU:
                    tmp.mult(TSUN[idst],TSPN[isrc]);
                    break;
                case UTOP:
                    tmp.mult(TSPN[idst],TSUN[isrc]);
                    break;
            }
            trace(tmpd,tmp);
            tmpd *= 2 ;  // because the normalization of SUN generators is 
                         // Tr(T^2) = 0.5
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
    
    auto dstAdjIndex = [dest_group_algebra_dim](int i, int j){
        return i*dest_group_algebra_dim + j;
    };
    auto srcAdjIndex = [source_group_algebra_dim](int i, int j){
        return i*source_group_algebra_dim + j;
    };

    stringstream ss;
    // src = spn, dst = sun:
    // "#define _spntnosun_adj(sunadjm,spnadjm) \\\n";
    // src = sun, dst = spn:
    // "#define _suntnospn_adj(spnadjm,sunadjm) \\\n";
    ss << "#define _" << src <<  "to" << dst << "_adj("<< dst <<"adjm,"<< src <<"adjm) \\\n";

    for(int isu1 = 0; isu1 < dest_group_algebra_dim; ++isu1)
        for(int isu2 = 0; isu2 < dest_group_algebra_dim; ++isu2){
            int rowIndex = dstAdjIndex(isu1,isu2);
            ss << "\t"<< dst << "adjm[" << rowIndex << "]=";
            map<int,double> newRow; // 
            auto row1 = dst_from_src.at(isu1);
            auto row2 = dst_from_src.at(isu2);
            for(auto e1 : row1) for(auto e2 : row2){
                int c1,c2;
                double w1,w2;
                tie(c1,w1) = e1;
                tie(c2,w2) = e2;
                int colIndex = srcAdjIndex(c1,c2);
                newRow[colIndex] += w1*w2; 
            }
            for(auto entry: newRow){
                int idx;
                double w;
                tie(idx,w) = entry;
                ss << "+(" << setprecision(15) << w << "*" << src << "adjm["
                    << idx << "])";
            }
        ss << ";\\\n";
        }
    ss << "\n\n";
    delete[] TSUN;
    delete[] TSPN;
    return ss.str();
}

string write_spntosun_adj(){
    return write_conversions_adj_tensors(PTOU);

}
string write_suntospn_adj(){
    return write_conversions_adj_tensors(UTOP);
}


void write_acf(string progname, string acfname, int N){ // algebra converter file
	group::N = N;
    ofstream acf(acfname.c_str()); 
    acf << "/** This file is produced automatically by " << progname << " */" << endl; 
    acf << "#define N " << N << endl << endl;
    acf << write_spntosun_adj();
    acf << write_suntospn_adj();
    acf << write_spntosun_algebra();
    acf << write_suntospn_algebra();
    acf.close();
}

void write_mtp(string progname, string acfname, int N, string mtpname){
    ofstream mtp(mtpname.c_str()); // macro test program

    mtp << "/** This file is produced automatically by " << progname << " */" << endl; 
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
    mtp << "              cout << \"Problem With component \" << setw(4) << i << \" \" ;" << endl;
    mtp << "              cout << spnAlgebraVectorCheck[i] << \" vs  \" << spnAlgebraVector[i];" << endl;
    mtp << "              cout << endl;" << endl;
    mtp << "              ok=false;" << endl;
    mtp << "          };" << endl;
    mtp << "      };" << endl;
    mtp << "      if(ok) cout << \"All seems ok for the algebra.\" << endl; " << endl;
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
    mtp << "      if(ok) cout << \"All seems ok for the adjmatrix.\" << endl; " << endl;
    mtp << "      else exit(1);" << endl;
    mtp << "  }" << endl;
    mtp << "  return 0; " << endl;
    mtp << "}" << endl;
    mtp.close();

}


