#include"nlib.h"
#include"math_nishi.h"
#include<Eigen/Dense>  // Eigen/Core, Geometry, LU, Cholesky, SVD, QR and Eigen values

using namespace Eigen;

float dihedral_4(Vector3f,Vector3f,Vector3f,Vector3f);

int dihedralnishi( Inp_nishi inp1 ){
  cout<<"Program dihedralnishi( inp1 ) starts!!!\n";
  cout<<"input file name = "<<inp1.filename <<endl;

// ############# READ PARAMETERS ########################################
/* (1) parameter settings
 * 
 */
  cout<<endl<<"--- READING PARAMETERS --- \n";
  string codname = inp1.read("COD1");
  string pdbname = inp1.read("REFPDBNAME");
  string selatom = inp1.read("SELATOM");
  int stride = atoi(inp1.read("STRIDE").c_str());
  int startframe = atoi( inp1.read("STARTFRAME").c_str() ) - 1 ;
  if(stride <= 0){
    return -1;
  }

// ############# LOAD COD FILES ########################################
/* (2) load trajectories by tra_nishi
 * 
 */
  cout<<"loading\n"<<codname<<" and "<<pdbname<<endl;
  tra_nishi* tra1;
  tra1 = new tra_nishi(codname.c_str(), pdbname.c_str(), stride, selatom);
  cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
  cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
  cout<<"TOTAL SELECTED ATOM = "<<tra1->total_sel<<endl;
  //unsigned int frame = tra1->total_step;

/* (2) setting of region
 *     from intra_start to intra_end (internal number of atom)
 */
  //cout<<endl<<"REPORT> (2) specify the region \n";
  cout<<endl<<"--- RESIDUE RANGE --- \n";
  string startchain = inp1.read("STARTCHAIN");
  //char endchain = inp1.read("ENDCHAIN").c_str()[0];
  int startres = atoi(inp1.read("STARTRES").c_str());
  //int endres = atoi(inp1.read("ENDRES").c_str());

// ############# CALCULATE phi psi ########################################
/* (x) calculate phi psi
 * 
 */

  vector<float> phi, psi; 

  //tra1->pdb1->write_pdb("zzz.pdb");
  //tra1->write_cod("zzz.pdb");
  int buf1,buf2,buf3,buf4,buf5;
  buf1 = search_sel( *tra1->pdb1, startchain, startres -1, "C", selatom);
  buf2 = search_sel( *tra1->pdb1, startchain, startres, "N", selatom);
  buf3 = search_sel( *tra1->pdb1, startchain, startres, "CA", selatom);
  buf4 = search_sel( *tra1->pdb1, startchain, startres, "C", selatom);
  buf5 = search_sel( *tra1->pdb1, startchain, startres +1, "N", selatom);
  
  Vector3f r1, r2, r3, r4, r5; // initial value is (1, 0, 0)
/*  for( unsigned int n = startframe; n < tra1->total_step; n++){

    //cout<<"search_sel( "<< startchain<<" "<<startres -1 <<" C"<<selatom<<"\n";
    //cout<<"buf1 = "<<buf1<<endl;
    r1 << tra1->cordx[n*tra1->total_sel+buf1], tra1->cordy[n*tra1->total_sel+buf1], tra1->cordz[n*tra1->total_sel+buf1];
    r2 << tra1->cordx[n*tra1->total_sel+buf2], tra1->cordy[n*tra1->total_sel+buf2], tra1->cordz[n*tra1->total_sel+buf2];
    r3 << tra1->cordx[n*tra1->total_sel+buf3], tra1->cordy[n*tra1->total_sel+buf3], tra1->cordz[n*tra1->total_sel+buf3];
    r4 << tra1->cordx[n*tra1->total_sel+buf4], tra1->cordy[n*tra1->total_sel+buf4], tra1->cordz[n*tra1->total_sel+buf4];
    r5 << tra1->cordx[n*tra1->total_sel+buf5], tra1->cordy[n*tra1->total_sel+buf5], tra1->cordz[n*tra1->total_sel+buf5];

    //float dih = dihedral_4(r1,r2,r3,r4);
    phi.push_back( dihedral_4(r1,r2,r3,r4) );
    psi.push_back( dihedral_4(r2,r3,r4,r5) );

    //cout<<"dih = "<<dih<<endl;

  }
*/
/*  calculate dihedral of reference pdb file and cout 
*/
/*  cout<<"calculate dihedral of reference pdb file and cout \n";
  cout<<"please put these dihedral angles into output-file, if you want to see them in the same time.\n";
  buf1 = search_sel( *tra1->pdb1, startchain, startres -1, "C", "all");
  buf2 = search_sel( *tra1->pdb1, startchain, startres, "N", "all");
  buf3 = search_sel( *tra1->pdb1, startchain, startres, "CA","all");
  buf4 = search_sel( *tra1->pdb1, startchain, startres, "C", "all");
  buf5 = search_sel( *tra1->pdb1, startchain, startres +1, "N", "all");
  r1 << tra1->pdb1->coox[buf1], tra1->pdb1->cooy[buf1], tra1->pdb1->cooz[buf1];
  r2 << tra1->pdb1->coox[buf2], tra1->pdb1->cooy[buf2], tra1->pdb1->cooz[buf2];
  r3 << tra1->pdb1->coox[buf3], tra1->pdb1->cooy[buf3], tra1->pdb1->cooz[buf3];
  r4 << tra1->pdb1->coox[buf4], tra1->pdb1->cooy[buf4], tra1->pdb1->cooz[buf4];
  r5 << tra1->pdb1->coox[buf5], tra1->pdb1->cooy[buf5], tra1->pdb1->cooz[buf5];
  cout<<"phi psi \n";
  cout<<dihedral_4(r1,r2,r3,r4)<<" "<<dihedral_4(r2,r3,r4,r5)<<endl;
*/
/*  output file
*/
/*  string outfile = inp1.read("OUTFILE").c_str();
  ofstream ofs;
  ofs.open( outfile.c_str() );
  int frame = tra1->total_step - startframe;
  for(int n=0;n<frame;n++){
    ofs<<phi[n]<<"   "<<psi[n]<<"   "<<n+1<<endl;
  }
  ofs.close();
  cout<<"output "<<outfile<<" (for graph) \n";
*/
// ############# CALCULATE phi psi moura ########################################
/* (x) moura  
 * 210-222 all ca 
 */
  int res1 = atoi(inp1.read("RES1").c_str());
  int res2 = atoi(inp1.read("RES2").c_str());
  int res3 = atoi(inp1.read("RES3").c_str());
  int res4 = atoi(inp1.read("RES4").c_str());
  vector<float> pseudo_dih;
  buf1 = search_sel( *tra1->pdb1, startchain, res1, "CA",selatom );
  buf2 = search_sel( *tra1->pdb1, startchain, res2, "CA",selatom );
  buf3 = search_sel( *tra1->pdb1, startchain, res3, "CA",selatom );
  buf4 = search_sel( *tra1->pdb1, startchain, res4, "CA",selatom   );
  for( unsigned int n = startframe; n < tra1->total_step; n++){
    r1 << tra1->cordx[n*tra1->total_sel+buf1], tra1->cordy[n*tra1->total_sel+buf1], tra1->cordz[n*tra1->total_sel+buf1];
    r2 << tra1->cordx[n*tra1->total_sel+buf2], tra1->cordy[n*tra1->total_sel+buf2], tra1->cordz[n*tra1->total_sel+buf2];
    r3 << tra1->cordx[n*tra1->total_sel+buf3], tra1->cordy[n*tra1->total_sel+buf3], tra1->cordz[n*tra1->total_sel+buf3];
    r4 << tra1->cordx[n*tra1->total_sel+buf4], tra1->cordy[n*tra1->total_sel+buf4], tra1->cordz[n*tra1->total_sel+buf4];
    pseudo_dih.push_back( dihedral_4(r1,r2,r3,r4) );
  }

  string outfile2 = inp1.read("OUTFILE2").c_str();
  ofstream ofs2;
  ofs2.open( outfile2.c_str() );
  for(int n=0;n<pseudo_dih.size();n++){
    ofs2<<pseudo_dih[n]<<endl;
    //cout<<pseudo_dih[n]<<endl;
  }
  ofs2.close();
  cout<<"output "<<outfile2<<" (for graph) \n";

  return 0;
}


float dihedral_4(Vector3f r1,Vector3f r2,Vector3f r3,Vector3f r4){
  Vector3f r12 = r2-r1, r23 = r3-r2, r34 = r4-r3;
  float sin = r12.dot(r23.cross(r34)) * r12.norm() / r12.cross(r23).norm() / r23.cross(r34).norm();
  float cos = r12.cross(r23).dot(r23.cross(r34)) / r12.cross(r23).norm() / r23.cross(r34).norm();
  float sign = sin / fabs(sin); //sign ( + or - ) of sin
  float t = sign * acos(cos); //torsion anglr t (tau) (rad)
  t = t * 180 / PI; //rad -> degree

  return t;
}


