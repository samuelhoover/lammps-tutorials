#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include<string>
//Polystyrene sulphonate is polyanion
using namespace std;

class vect{
	public:
		int r[3];
		void assign(int *);
};

class co{
	public:
		int index;
		int molec, type;
		int charge;
		int r[3];

		void assign(int , int, int, int, int *);
		void Translate(int *);
		void get_value(int *);
		void print();
};

class bo{
	public:
		int index, type;
		int a1, a2;
		void assign(int, int, int, int);
		void print();
};

class ang{
	public:
		int index, type;
		int a1, a2, a3;
		void assign(int, int, int, int,int);
		void print();
};

co *chain1, *cation, *anion;
vect *transf;
bo *bond1;
ang *angle1;
int Nbond, Nchain, Nangle, Ncation, Nanion, Npolym;
double QFAC, salt_conc, polym_conc;
char param_file[200];
int box_max[3], box_min[3];
int max_lat[3];

/******co functions *****************/

void vect::assign(int *r1)
{
	for(int i=0;i<3;i++) r[i]=r1[i]*max_lat[i];
}

void co::assign(int ind, int m, int t, int c, int *r1)
{
	int i;
	index=ind;
	molec=m;
	type=t;
	charge=c;
	for(i=0;i<3;i++)
		r[i]=r1[i];
}


void co::Translate(int *vec)
{
	int i;
	for(i=0;i<3;i++)
		r[i]+=vec[i];
}

void co::get_value(int *r1)
{
	int i;
	for(i=0;i<3;i++)
		r1[i]=r[i];
}

void co::print()
{
	int i;
	char c;
	c=type+64;
//	out<<index<<" "<<molec<<" "<<type<<" "<<charge<<" "; //conversion factor for charge
//	cout<<c<<" ";
	for(i=0;i<3;i++){
//		out<<r[i]<<" ";
//		cout<<r[i]<<" ";
		}
//	out<<endl;
//	cout<<endl;
}
/*************bond ***********/
void bo::assign(int b1, int b2, int a11, int a21)
{
	index=b1; type=b2;
	a1=a11; a2=a21;	
}

void bo::print()
{
//	out<<index<<" "<<type<<" "<<a1<<" "<<a2<<endl;
}

/**************angle **************/
void ang::assign(int b1, int b2, int c1, int c2, int c3)
{
	index=b1;
	type=b2;
	a1=c1; a2=c2; a3=c3;
}

void ang::print()
{
//	out<<index<<" "<<type<<" "<<a1<<" "<<a2<<" "<<a3<<endl;
}

int get_dist_2(int *r1, int *r2)
{
	int i, dist=0;
	for(i=0;i<3;i++)dist+=(r1[i]-r2[i])*(r1[i]-r2[i]);
	return dist;
}

int is_occupied(int index, int *r)
{
	int i, dist_2;

	//if previous coordinate is occupied
	for(i=0;i<index;i++)
	{
		dist_2=get_dist_2(r, chain1[i].r);
		if(dist_2==0) return 1;
	}

	return 0;

}

int is_overlap_chains(int *r1)
{
	int k, l,j, r2[3], dist_2;
	//check for all monomers
	for(k=0;k<Nchain;k++){
		for(l=0;l<2*Npolym;l++){
			for(j=0;j<3;j++)
				r2[j]=chain1[k].r[j]+transf[l].r[j];
			 dist_2=get_dist_2(r1,r2);
			if(dist_2==0) return 1;
		}
	}
	return 0;		
}

int is_overlap_ions(int ind1, int ind2, int *r1)
{
	int k, dist_2;
	for(k=0;k<ind1;k++){
		dist_2=get_dist_2(r1, cation[k].r);
		if(dist_2==0) return 1;
	}
	for(k=0;k<ind2;k++){
		dist_2=get_dist_2(r1, anion[k].r);
		if(dist_2==0) return 1;
	}
	return 0;
}

int get_coord_ion(int index, int tp, int *r)
{
	int i;
	int r1[3], isoverlap=1, ind1;
	int cut_off;

	cut_off=Ncation/Npolym+1;

	while(isoverlap){
	 //cation
		ind1=index/cut_off;
	if(tp) //anion
		ind1+=Npolym;

	for(i=0;i<3;i++)
		r1[i]=rand()%(2*max_lat[i]+1)-max_lat[i]+transf[ind1].r[i];
	//check for overlap with monomers
	isoverlap=is_overlap_chains(r1);

	//check for the previous cations
	if(!isoverlap){
	if(tp==0){  //if cation
	if(index>0)
		isoverlap=is_overlap_ions(index, 0, r1);
		}
	else {  //if anion
		isoverlap=is_overlap_ions(Ncation,index,r1);	
	}
	}
	}

	for(i=0;i<3;i++) r[i]=r1[i];
	return 1;
}

void print_chain(int tm)
{
	int i,j;
	ofstream out;
	out.open("dump.lammpstrj");

	out<<"ITEM: TIMESTEP"<<endl;
	out<<tm<<endl;
	out<<"ITEM: NUMBER OF ATOMS\n";
	out<<Nchain<<endl;
	out<<"ITEM: BOX BOUNDS\n";

	for(i=0;i<3;i++)
		out<<box_min[i]<<"\t"<<box_max[i]<<endl;
	out<<"ITEM: ATOMS id type xu yu zu\n";

	for(i=0;i<Nchain;i++)
	{
		out<<chain1[i].index<<"\t"<<chain1[i].type<<"\t";
			for(j=0;j<3;j++)
			out<<chain1[i].r[j]<<"\t";
		out<<endl;
	}

	out.close();
}

/**************create single chain ***************/
int create_single_chain()
{
	int i,j,k,l=1;
	int mol=1,tp=1, index=1;
        int ch=-1;	
	int r[3]={0,0,Nchain/2}, r1[3]={0,0,0}, tmp;
	int isoccupied;
	for(i=0;i<3;i++)
	max_lat[i]=0;
	

	
	chain1=new co[Nchain];
	for(i=0;i<Nchain;i++){
		chain1[i].assign(i+1, mol, tp, ch, r);
		r[2]-=1;
	}

	for(i=0;i<Nbond;i++){
	bond1[i].assign(i+1,1, i+1, i+2);
	}

	for(i=0;i<Nangle;i++){
	angle1[i].assign(i+1,1, i+1, i+2, i+3);
	}

	for(i=0;i<3;i++)
		max_lat[i]=Nchain/8;

	print_chain(0);
	cout<<"exiting create single chain\n";
	return 0;

}

int get_trans_vect(int index1, int *trans)
{
	int i=0,j=0,k=0, l=0;
	int lim[3];
	double tmp;

	if(index1==0) {
	for(j=0;j<3;j++) trans[j]=0;
	return 1;
	}

	for(i=0;i<3;i++){
		tmp=(box_max[i]-box_min[i]-2*max_lat[i])/(2*max_lat[i]);
		lim[i]=ceil(tmp);
		cout<<"limit = "<<lim[i]<<endl;
		}
	
	for(i=-lim[0]/2;i<lim[0]/2;i++){
		for(j=-lim[1]/2;j<lim[1]/2;j++){
			for(k=-lim[2]/2;k<lim[2]/2;k++){
			l++;
			if(l==index1){
			trans[0]=i; trans[1]=j; trans[2]=k;
			return 1;
			}
			}
			}
			}
			return 0;

}



int create_data_file()
{
	ofstream out;
	int i,j,k,l, trans[3], molec_index, atom_index, tp=2, ch=1;
	Nangle=0;
	out.open("data.complexation");
	
	out<<"LAMMPS Description\n\n";
	out<<Npolym*Nchain<<" atoms\n";
	out<<Npolym*Nbond<<" bonds\n";
	out<<Npolym*Nangle<<" angles\n";
	out<<"0 dihedrals\n0 impropers\n\n";
	out<<"1 atom types\n1 bond types\n0 angle types\n0 dihedral types\n0 improper types\n\n";

	out<<box_min[0]<<" "<<box_max[0]<<" xlo xhi\n";
	out<<box_min[1]<<" "<<box_max[1]<<" ylo yhi\n";
	out<<box_min[2]<<" "<<box_max[2]<<" zlo zhi\n";
	
	out<<"\nMasses\n\n1 1.0\n\n";

	out<<"Atoms\n\n";
	
	//Create the polymer chains by translating orginal chain
	transf=new vect[Npolym];
	//Create Polyanion first
	for(k=0;k<Npolym;k++)
	{
	l=get_trans_vect( k, trans);
	transf[k].assign(trans);
	for(i=0;i<Nchain;i++)
	{
		out<<chain1[i].index+k*Nchain<<" "<<chain1[i].molec+k<<" "<<chain1[i].type<<" "<<chain1[i].charge*QFAC<<" ";
		for(j=0;j<3;j++){
			out<<chain1[i].r[j]+trans[j]*max_lat[j]<<" ";
		}
		out<<endl;
	}
	}


	molec_index=Npolym+1; atom_index=Npolym*Nchain+1;
	//create the cations
	ch=1;  //cation charge
	trans[0]=1; trans[1]=0; trans[2]=Nchain/2;
/*	for(i=0;i<Ncation;i++)
	{
	
	cation[i].assign(atom_index++, molec_index++, tp, ch, trans); 
	out<<cation[i].index<<" "<<cation[i].molec<<" "<<cation[i].type<<" "<<cation[i].charge*QFAC<<" ";
	for(j=0;j<3;j++)
		out<<cation[i].r[j]<<" ";
	out<<endl;	
	
	if(i%2==0)
	{
		trans[0]=0; trans[1]=1;
	}
	else{
		trans[0]=1; trans[1]=0;
	}
	trans[2]-=1;

	}
*/


	out<<"\nBonds\n\n";

	for(k=0;k<Npolym;k++)
	{
	for(i=0;i<Nbond;i++)
	{
		out<<bond1[i].index+k*Nbond<<" "<<bond1[i].type<<" "<<bond1[i].a1+k*Nchain<<" "<<bond1[i].a2+k*Nchain<<endl;	
	}
	}



	out.close();
	return 0;		
}

int create_in_file(int rnd1, int res)
{	
	ofstream out;
	ifstream in;
	double tmp[5];
	char str[200];
	int i,j,k;

	out.open("in.complexation");
	in.open(param_file);

	out<<"###################### polyelectrolyte complexation ######################\n";
	out<<"clear\n\n#Initialization\n#------------------------------------\n";
	out<<"units\t\t\t\tlj\ndimension\t\t\t3\natom_style\t\t\tfull\nboundary\t\t\tp p p\n\n";
	out<<"#Atom Definition\n";
	out<<"#-------------------------------------\n";
	if(!res)
	out<<"read_data\t\t\tdata.complexation\n\n";
	else if(res==1)
	out<<"read_restart\t\t\trestart.equil\n\n";
	else
	out<<"read_restart\t\t\trestart.rerun."<<res-1<<"\n\n";


	out<<"#Force definition\n#----------------------------------------------\n#Bond definition\n";
	out<<"bond_style\t\t\tharmonic\n";
	in>>str;
	for(i=0;i<3;i++) in>>tmp[i];
	out<<"bond_coeff\t\t\t";
	for(i=0;i<3;i++)
		out<<tmp[i]<<" ";
	out<<"\n#Pair definition\n";
	out<<"special_bonds\tlj 0 1 1 \n";
	in>>str>>tmp[0];
	if(res)
	out<<"dielectric\t\t"<<tmp[0]<<"\n";
	in>>str>>j;
	in>>str>>tmp[0]>>tmp[1];

	out<<"pair_style\t\t lj/cut "<<tmp[0]<<"\n";

	out<<"pair_modify \t\t shift yes\n";
	for(i=0;i<j;i++){
	in>>str;
	out<<"pair_coeff\t";
	for(k=0;k<5;k++){
		in>>tmp[k];
		out<<tmp[k]<<" ";
		}
		out<<"\n";

	}
	//Read the angle	
	in>>str;
	for(i=0;i<3;i++) in>>tmp[i];
//	out<<"angle_style cosine/squared\n";
//	out<<"angle_coeff\t\t\t";
//	for(i=0;i<3;i++)
//	out<<tmp[i]<<"\t";
//	out<<endl;
	in>>str>>tmp[0];
	cout<<"PPM "<<tmp[0]<<endl;
//	if(res)
//	out<<"kspace_style\t\tpppm\t\t"<<tmp[0]<<endl;

	out<<"\n#Timestep etc\n#--------------------------------------------\n";
	in>>str>>tmp[0]>>str>>tmp[1]>>str>>tmp[2];
	if(!res) tmp[0]/=100.0;
	out<<"timestep\t\t"<<tmp[0]<<"\nrun_style\t\tverlet\nvelocity\t\tall\tcreate\t"<<tmp[1]<<"\t"<<rnd1<<endl;
	out<<"\n#Fix\n#---------------------------------------\n";

	if(res==0 || res==1)
	out<<"group polymers type 1\n";

	out<<"fix 1 all nve\n";
	out<<"fix 2 all langevin "<<tmp[1]<<"\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t"<<rnd1<<endl;

	if(res==0 )
	out<<"fix 3 polymers recenter INIT INIT INIT units box\n";

	out<<"\n#Dump\n#------------------------------------------------\n";
	out<<"thermo_style\t\t custom step temp press pe evdwl ecoul ebond ke etotal enthalpy\n";
	in>>str>>tmp[3]>>str>>tmp[4]>>str>>tmp[0];
	out<<"thermo \t"<<tmp[3]<<endl;
	out<<"dump\t3 all custom "<<tmp[3]<<" dump.lammpstrj "<<"id mol type xu yu zu\n";
	//out<<"dump\t3 all custom "<<tmp[3]<<" traj.lammpstrj "<<"id mol type x y z\n";
	
	out<<"run\t\t" << fixed << int(tmp[4]) << endl;

		if(!res)
	out<<"write_restart restart.equil\n";
	else
	out<<"write_restart restart.rerun."<<res<<"\n"; 

	out<<"#--------------End of Input file ------------------\n";
	
	QFAC=0;
	out.close();
	in.close();
	
	return 0;
}

void cleanup()
{
	if(chain1!=NULL)  delete []chain1;
	if(cation!=NULL)  delete []cation;
	if(anion!=NULL)  delete []anion;
	if(bond1!=NULL)  delete []bond1;
	if(angle1!=NULL) delete []angle1;
	if(transf!=NULL) delete []transf;
}

void read_config(char *str)
{
	int i;
	double tmp;
	char s1[500], s2[500];
	ifstream in;
	in.open(str);
	
	for(i=0;i<3;i++)
		in.getline(s1,500);
	in>>s1>>param_file;
	cout<<param_file<<endl;

	for(i=0;i<5;i++)
	in.getline(s1,500);

	//Make sure the box limits are integers
	for(i=0;i<3;i++){
		in>>s1>>s2;
		box_min[i]=atof(s1); box_max[i]=atof(s2);
	}

	for(i=0;i<3;i++)
	in.getline(s1,500);

	//Initialize the chain and bond parameters 
	in>>s1>>s1;
	Nchain=atoi(s1);
	cout<<"NChain = "<<Nchain<<endl;
	Nbond=Nchain-1;
	bond1=new bo[Nbond];
	Nangle=Nchain-2;
	angle1=new ang[Nangle];
	in>>s1>>s1;
	polym_conc=atof(s1);

	//Reset the box dimension
/*	for(i=0;i<3;i++)
	{
		if(i==0 || i==1){
		box_min[i]= -Nchain/4; box_max[i]= Nchain/4;}

		else{
			box_min[i]=-Nchain/2-10;
			box_max[i]=Nchain/2+10;}

	}

*/
	//Number of polymer chain.. convert from mole/L to molecule/nm^3 then multiply by box in nm^3
		tmp=0.6023*polym_conc;
		for(i=0;i<3;i++)
			tmp*=(box_max[i]-box_min[i])*0.25; //0.25 nm is the distance between consecutive charges on the polymer
		Npolym	=ceil(tmp);
		if(Npolym<1) Npolym=1;
	
		cout<<"N_polym= "<<Npolym<<endl;

	for(i=0;i<3;i++)
	in.getline(s1,500);
	in>>s1>>s1;
	salt_conc=atof(s1);
	cout<<"salt conc= "<<s1<<endl;

	Ncation=Nchain*Npolym; Nanion=0;
	//Add salt concentration by converting real and lj units
	tmp=0.6023*salt_conc;
	for(i=0;i<3;i++)
		tmp*=(box_max[i]-box_min[i])*0.25;
	int Nsalt;
	Nsalt=ceil(tmp);
	cout<<"N_salt= "<<Nsalt<<endl;
	Ncation+=Nsalt; Nanion+=Nsalt;
	cout<<"exiting read config\n";
	cation=new co[Ncation];
	anion=new co[Nanion];
		
	in.close();
}

/*************main *********************/
//./a.out config.in <random number> <equil=0 otherwise 1>
int main(int argc, char **argv)
{
	int i,j;
	//reading the configuration file
	if(argc>3){
	read_config(argv[1]);
	i=atoi(argv[2]); //random number as an input for n number of ensemble
	j=atoi(argv[3]); //For equilibriation step set 0 otherwise 1
	}
	else{
		cout<<"Enter config file name random number and step number\n";
		return 0;
	}
	srand(i); //Random number generation
	create_single_chain(); //Create 1 chain random walk with excluded volume effect 
	create_in_file(i,j); //Create in file with j value used for equilibration or electrostatics
	if(!j)
	create_data_file(); //create data file only for equilibration else read restart from equilibrated run
	cleanup();
	return 0;
}

