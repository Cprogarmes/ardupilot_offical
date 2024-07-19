
/* PSINS(Precise Strapdown Inertial Navigation System) C++ algorithm source file PSINS.cpp

Copyright(c) 2015-2023, by YanGongmin, All rights reserved.
Northwestern Polytechnical University, Xi'an, P.R.China.
Date: 17/02/2015, 19/07/2017, 11/12/2018, 27/12/2019, 12/12/2020, 22/11/2021, 17/10/2022, 23/08/2023
*/

#include "PSINS.h"

const CVect3 O31(0.0), One31(1.0), I31Z(0,0,1.0), Ipos(1.0/RE,1.0/RE,1.0), posNWPU=LLH(34.246048,108.909664,380); //NWPU-Lab-pos
const CQuat  qI(1.0,0.0,0.0,0.0);
const CMat3  I33(1,0,0, 0,1,0, 0,0,1), O33(0,0,0, 0,0,0, 0,0,0), One33(1.0);
const CVect  On1(MMD,0.0), O1n=~On1, Onen1(MMD,1.0);
CVect3       Vrbs(0.0);
CGLV		 glv;
int			 psinslasterror = 0;
int			 psinsstack0 = 0, psinsstacksize = 0;

//***************************  class CGLV  *********************************/
CGLV::CGLV(double Re0, double f0, double g00, double wie0)
{
	this->Re = Re0; this->f = f0; this->g0 = g00; this->wie = wie0;
	Rp = (1-f0)*Re0;
	e = sqrt(2*f0-f0*f0); e2 = e*e;
	ep = sqrt(Re0*Re0-Rp*Rp)/Rp; ep2 = ep*ep;
    mg = g00/1000.0;
    ug = mg/1000.0;
    deg = PI/180.0;
    min = deg/60.0;
    sec = min/60.0;
    ppm = 1.0e-6;
    hur = 3600.0;
	dps = deg/1.0;
    dph = deg/hur;
    dpsh = deg/sqrt(hur);
    dphpsh = dph/sqrt(hur);
    dph2 = dph/hur;
	dphpg = dph/g00;
    ugpsHz = ug/sqrt(1.0);
    ugpsh = ug/sqrt(hur);
	ugpg2 = ug/g00/g00;
    mpsh = 1/sqrt(hur); 
    mpspsh = 1/1/sqrt(hur);
    ppmpsh = ppm/sqrt(hur);
    secpsh = sec/sqrt(hur);
#ifdef PSINS_IO_FILE
	t0 = clock();
#endif
}

//***************************  class CVect3  *********************************/
CVect3::CVect3(void)
{
}

CVect3::CVect3(double xyz)
{
	i=j=k=xyz;
}

CVect3::CVect3(double xx, double yy, double zz)
{
	i=xx, j=yy, k=zz;
}

CVect3::CVect3(const double *pdata)
{
	i=*pdata++, j=*pdata++, k=*pdata;
}

CVect3::CVect3(const float *pdata)
{
	i=*pdata++, j=*pdata++, k=*pdata;
}

BOOL IsZero(const CVect3 &v, double eps)
{
	return (v.i<eps&&v.i>-eps && v.j<eps&&v.j>-eps && v.k<eps&&v.k>-eps);
}

BOOL IsZeroXY(const CVect3 &v, double eps)
{
	return (v.i<eps&&v.i>-eps && v.j<eps&&v.j>-eps);
}

BOOL IsNaN(const CVect3 &v)
{
	return 0; //(_isnan(i) || _isnan(j) || _isnan(k));
}

CVect3 CVect3::operator+(const CVect3 &v) const
{
	return CVect3(this->i+v.i, this->j+v.j, this->k+v.k);
}

CVect3 CVect3::operator-(const CVect3 &v) const
{
	return CVect3(this->i-v.i, this->j-v.j, this->k-v.k);
}

CVect3 CVect3::operator*(const CVect3 &v) const
{
	//stacksize();
	return CVect3(this->j*v.k-this->k*v.j, this->k*v.i-this->i*v.k, this->i*v.j-this->j*v.i);
}
	
CVect3 CVect3::operator*(double f) const
{
	return CVect3(i*f, j*f, k*f);
}

CVect3 CVect3::operator*(const CMat3 &m) const
{
	return CVect3(i*m.e00+j*m.e10+k*m.e20,i*m.e01+j*m.e11+k*m.e21,i*m.e02+j*m.e12+k*m.e22);
}
	
CVect3 CVect3::operator/(double f) const
{
	return CVect3(i/f, j/f, k/f);
}

CVect3 CVect3::operator/(const CVect3 &v) const
{
	return CVect3(i/v.i, j/v.j, k/v.k);
}

CVect3& CVect3::operator=(double f)
{
	i = j = k = f;
	return *this;
}

CVect3& CVect3::operator=(const double *pf)
{
	i = *pf++, j = *pf++, k = *pf;
	return *this;
}

CVect3& CVect3::operator+=(const CVect3 &v)
{ 
	i += v.i, j += v.j, k += v.k;
	return *this;
}

CVect3& CVect3::operator-=(const CVect3 &v)
{ 
	i -= v.i, j -= v.j, k -= v.k;
	return *this;
}

CVect3& CVect3::operator*=(double f)
{ 
	i *= f, j *= f, k *= f;
	return *this;
}

CVect3& CVect3::operator/=(double f)
{ 
	i /= f, j /= f, k /= f;
	return *this;
}

CVect3& CVect3::operator/=(const CVect3 &v)
{ 
	i /= v.i, j /= v.j, k /= v.k;
	return *this;
}

CVect3 operator*(double f, const CVect3 &v)
{
	return CVect3(v.i*f, v.j*f, v.k*f);
}
	
CVect3 operator-(const CVect3 &v)
{
	return CVect3(-v.i, -v.j, -v.k);
}

double& CVect3::operator()(int r)
{
	return (&this->i)[r];
}

double crossXY(const CVect3 &v1, const CVect3 &v2)
{
	return v1.i*v2.j-v1.j*v2.i;
}

CMat3 vxv(const CVect3 &v1, const CVect3 &v2)
{
	return CMat3(v1.i*v2.i, v1.i*v2.j, v1.i*v2.k, 
				 v1.j*v2.i, v1.j*v2.j, v1.j*v2.k, 
				 v1.k*v2.i, v1.k*v2.j, v1.k*v2.k);
}

CVect3 sqrt(const CVect3 &v)
{
	return CVect3(sqrt(v.i), sqrt(v.j), sqrt(v.k));
}

CVect3 pow(const CVect3 &v, int k)
{
	CVect3 pp = v;
	for(int i=1; i<k; i++)
	{
		pp.i *= v.i, pp.j *= v.j, pp.k *= v.k;
	}
	return pp;
}

CVect3 abs(const CVect3 &v)
{
	CVect3 res;
	res.i = v.i>0.0 ? v.i : -v.i;
	res.j = v.j>0.0 ? v.j : -v.j;
	res.k = v.k>0.0 ? v.k : -v.k;
	return res;
}

CVect3 maxabs(const CVect3 &v1, const CVect3 &v2)
{
	CVect3 res1, res2;
	res1.i = v1.i>0.0 ? v1.i : -v1.i;  res2.i = v2.i>0.0 ? v2.i : -v2.i;  if(res1.i<res2.i) res1.i=res2.i;
	res1.j = v1.j>0.0 ? v1.j : -v1.j;  res2.j = v2.j>0.0 ? v2.j : -v2.j;  if(res1.j<res2.j) res1.j=res2.j;
	res1.k = v1.k>0.0 ? v1.k : -v1.k;  res2.k = v2.k>0.0 ? v2.k : -v2.k;  if(res1.k<res2.k) res1.k=res2.k;
	return res1;
}

double norm(const CVect3 &v)
{
	return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

double normlize(CVect3 *v)
{
	double n=norm(*v);
	v->i /= n, v->j /= n, v->k /= n;
	return n;
}

double normInf(const CVect3 &v)
{
	double i = v.i>0 ? v.i : -v.i,
		   j = v.j>0 ? v.j : -v.j,
		   k = v.k>0 ? v.k : -v.k;
	if(i>j)	return i>k ? i : k;
	else    return j>k ? j : k;
}

double normXY(const CVect3 &v)
{
	return sqrt(v.i*v.i + v.j*v.j);
}

double normXYInf(const CVect3 &v)
{
	double i = v.i>0 ? v.i : -v.i,
		   j = v.j>0 ? v.j : -v.j;
	return i>j ? i : j;
}

double dot(const CVect3 &v1, const CVect3 &v2)
{
	return (v1.i*v2.i + v1.j*v2.j + v1.k*v2.k);
}

CVect3 dotmul(const CVect3 &v1, const CVect3 &v2)
{
	return CVect3(v1.i*v2.i, v1.j*v2.j, v1.k*v2.k);
}

CVect3 dotdiv(const CVect3 &v1, const CVect3 &v2)
{
	return CVect3(v1.i/v2.i, v1.j/v2.j, v1.k/v2.k);
}

CQuat rv2q(const CVect3 &rv)
{
#define rvF1	(     2 * 1)		// define: Fk=2^k*k! 
#define rvF2	(rvF1*2 * 2)
#define rvF3	(rvF2*2 * 3)
#define rvF4	(rvF3*2 * 4)
#define rvF5	(rvF4*2 * 5)
	double n2 = rv.i*rv.i+rv.j*rv.j+rv.k*rv.k, c, f;
	if(n2<(PI/180.0*PI/180.0))	// 0.017^2 
	{
		double n4=n2*n2;
		c = 1.0 - n2*(1.0/rvF2) + n4*(1.0/rvF4);
		f = 0.5 - n2*(1.0/rvF3) + n4*(1.0/rvF5);
	}
	else
	{
		double n_2 = sqrt(n2)/2.0;
		c = cos(n_2);
		f = sin(n_2)/n_2*0.5;
	}
	return CQuat(c, f*rv.i, f*rv.j, f*rv.k);
}

CMat3 rv2m(const CVect3 &rv)
{
	return q2mat(rv2q(rv));
}

CMat3 askew(const CVect3 &v)
{
	return CMat3(0,  -v.k, v.j, 
				 v.k, 0.0,  -v.i,
				-v.j, v.i, 0);
}

CMat3 pos2Cen(const CVect3 &pos)
{
	double si = sin(pos.i), ci = cos(pos.i), sj = sin(pos.j), cj = cos(pos.j);
	return CMat3(	-sj, -si*cj,  ci*cj,  
					 cj, -si*sj,  ci*sj,  
					 0,   ci,     si      );	//Cen
}

CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts, CEarth *pEth)
{
	double sl, cl, sl2, sq, sq2, RMh, RNh, clRNh;
	if(pEth)
	{
		RMh = pEth->RMh; clRNh = pEth->clRNh;
	}
	else
	{
		sl=sin(pos0.i); cl=cos(pos0.i); sl2=sl*sl;
		sq = 1-glv.e2*sl2; sq2 = sqrt(sq);
		RMh = glv.Re*(1-glv.e2)/sq/sq2+pos0.k;
		RNh = glv.Re/sq2+pos0.k;    clRNh = cl*RNh;
	}
    CVect3 vn = pos1 - pos0;
    return CVect3(vn.j*clRNh/ts, vn.i*RMh/ts, vn.k/ts);
}

CVect3 pp2att(const CVect3 &pos1, const CVect3 &pos0)
{
	return vn2att(pp2vn(pos1, pos0, 1.0));
}

double sinAng(const CVect3 &v1, const CVect3 &v2, const CVect3 &v0)
{
	if(&v0!=&O31)
		return sinAng(v1-v0,v2-v0);
	if(IsZero(v1)||IsZero(v2)) return 0.0;
	return norm(v1*v2)/(norm(v1)*norm(v2));
}

double MagYaw(const CVect3 &mag, const CVect3 &att, double declination)
{
	CVect3 attH(att.i, att.j, 0.0);
	CVect3 magH = a2mat(attH)*mag;
	double yaw = 0.0;
//	if(attH.i<(80.0*DEG)&&attH.i>-(80.0*DEG))
	{
		yaw = atan2Ex(magH.i, magH.j) + declination;
		if(yaw>PI)       yaw -= _2PI;
		else if(yaw<-PI) yaw += _2PI;
	}
	return yaw;
}

CVect3 xyz2blh(const CVect3 &xyz)
{
	double s = normXY(xyz), theta = atan2(xyz.k*glv.Re, s*glv.Rp),
		s3 = sin(theta), c3 = cos(theta); s3 = s3*s3*s3, c3 = c3*c3*c3;
	if(s<(6378137.0*1.0*DEG))  return O31;
	double L = atan2(xyz.j, xyz.i), B = atan2(xyz.k+glv.ep2*glv.Rp*s3, s-glv.e2*glv.Re*c3),
		sB = sin(B), cB = cos(B), N = glv.Re/sqrt(1-glv.e2*sB*sB);
	return CVect3(B, L, s/cB-N);
}

CVect3 blh2xyz(const CVect3 &blh)
{
	double sB = sin(blh.i), cB = cos(blh.i), sL = sin(blh.j), cL = cos(blh.j),
		N = glv.Re/sqrt(1-glv.e2*sB*sB);
	return CVect3((N+blh.k)*cB*cL, (N+blh.k)*cB*sL, (N*(1-glv.e2)+blh.k)*sB);
}

CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos)
{
	return Vxyz*pos2Cen(pos);
}

CVect3 v2double(double f)
{
	CVect3 v;
	v.i = (int)(f+EPS);	f = (f-v.i)*1.0e8;  if(v.i>50000000.0) v.i -= 100000000.0;  // -45000000 ~ +45000000
	v.j = (int)(f+EPS);						if(v.j>50000000.0) v.j -= 100000000.0;
	v.k = 0.0;
	return v;
}

CVect3 v3double(double f)
{
	CVect3 v;
	v.i = (int)(f+EPS);	f = (f-v.i)*1.0e5;  if(v.i>50000.0) v.i -= 100000.0;  // -45000 ~ +45000
	v.j = (int)(f+EPS);	f = (f-v.j)*1.0e5;  if(v.j>50000.0) v.j -= 100000.0;
	v.k = (int)(f+EPS);						if(v.k>50000.0) v.k -= 100000.0;
	return v;
}

void v5double(double f, CVect3 &v1, CVect3 &v2)
{
	psinsassert(f>0.0);
	v1.i = (int)(f+EPS);	f = (f-v1.i)*1.0e3;   // 0 ~ +990
	v1.j = (int)(f+EPS);	f = (f-v1.j)*1.0e3;
	v1.k = (int)(f+EPS);	f = (f-v1.k)*1.0e3;
	v2.i = (int)(f+EPS);	f = (f-v2.i)*1.0e3;
	v2.j = (int)(f+EPS);
	v2.k = 0.0;
}

CVect3 sort(const CVect3 &v)
{
	CVect3 vtmp=v;
	if(vtmp.i<vtmp.j) swapt(vtmp.i,vtmp.j,double);
	if(vtmp.i<vtmp.k) swapt(vtmp.i,vtmp.k,double);
	if(vtmp.j<vtmp.k) swapt(vtmp.j,vtmp.k,double);
	return vtmp;
}

double median(const double &f1, const double &f2, const double &f3)
{
	return median(CVect3(f1,f2,f3));
}

double median(const CVect3 &v)
{
	double res;
	if(v.i>v.j) {
		if(v.j>v.k) res=v.j;
		else if(v.k>v.i) res=v.i;
		else res=v.k;
	}
	else { // v.i<=v.j
		if(v.j<v.k) res=v.j;
		else if(v.i>v.k) res=v.i;
		else res=v.k;
	}
	return res;
}

CVect3 dm2r(const CVect3 &v, int typ)
{
	if(typ==2) return CVect3(dm2r(v.i), dm2r(v.j), v.k);
	return CVect3(dm2r(v.i), dm2r(v.j), dm2r(v.k));
}

CVect3 v3mmm(const CVect3 &p1, const CVect3 &p2, const CVect3 &p3, CVect3 *pmm, BOOL isMean, const CVect3 &th)
{
	CVect3 res, mm;
	double *pf0=&res.i, *pm=&mm.i; const double *pf1=&p1.i, *pf2=&p2.i, *pf3=&p3.i;
	for(int i=0; i<3; i++,pf0++,pm++,pf1++,pf2++,pf3++) {
		if(*pf1>*pf2) {
			if(*pf2>*pf3) { *pf0=*pf2, *pm=*pf1-*pf3; }
			else if(*pf3>*pf1) { *pf0=*pf1, *pm=*pf3-*pf1; }
			else { *pf0=*pf3, *pm=*pf1-*pf2; }
		}
		else { // pf1<pf2
			if(*pf2<*pf3) { *pf0=*pf2, *pm=*pf3-*pf1; }
			else if(*pf1>*pf3) { *pf0=*pf1, *pm=*pf2-*pf3; }
			else { *pf0=*pf3, *pm=*pf2-*pf1; } 
		}
	}
	if(pmm) *pmm=mm;  // max-min
	if(isMean) {
		if(mm.i<th.i) res.i=(p1.i+p2.i+p3.i)/3.0;
		if(mm.j<th.j) res.j=(p1.j+p2.j+p3.j)/3.0;
		if(mm.k<th.k) res.k=(p1.k+p2.k+p3.k)/3.0;
	}
	return res;
}

CVect3 fopp(const CVect3 &a, const CVect3 &b, const CVect3 &c)  // point c to line a-b
{
	CVect3 ab=b-a;
	double u=dot(c-a,ab), n2=ab.i*ab.i+ab.j*ab.j+ab.k*ab.k;
	return n2<1.0e-20 ? a : a+ab*u/n2;
}

CVect3 fopp(const CVect3 &a, const CVect3 &b, const CVect3 &c, const CVect3 &d)  // point d to plane a-b-c, https://www.docin.com/p-2324716540.html
{
	CVect3 ABC=-inv(CMat3(a,b,c))*One31, ab=b-a, ac=c-a;
	return inv(CMat3(ab,ac,ABC))*CVect3(dot(d,ab),dot(d,ac),-1);
}

CVect3 tp2att(const CVect3 &a, const CVect3 &b, const CVect3 &c)  // a-b-c in couter-clockwise
{
	CVect3 d=fopp(a, b, c);
	CVect3 att = dv2att(c-d,a-d, CVect3(1,0,0),CVect3(0,1,0));  // Att of c-d-a(x-o-y)
	CVect3 vo = fopp(a, b, c, O31);  // orthocenter
//	CVect3 vo = (a+b+c)/3.0;  // center of gravity
	double afa=asin(sinAng(vo, b, a));
	return q2att(a2qua(att)*a2qua(CVect3(0,0,afa)));  // Att of *-fop-a(x-o-y)
}

CVect3 attract(const CVect3 &v, const CVect3 &th, const CVect3 &center)
{
	CVect3 dv=v-center;
	if(dv.i>th.i || dv.i<-th.i || dv.j>th.j || dv.j<-th.j || dv.k>th.k || dv.k<-th.k) return v;
	return CVect3(attract(v.i,th.i,center.i), attract(v.j,th.j,center.j), attract(v.k,th.k,center.k));
}

CVect3 ff2muxy(const CVect3 &f0, const CVect3 &f1, const char *dir0, const char *dir1)
{
    CVect3 v0 = f0, v1 = f1;
    if(dir0) IMURFU(&v0, 1, dir0);  
	if(dir1) IMURFU(&v1, 1, dir1);
    double n = norm(v0), f = cos(v0.j/n);  if(f<0.1) f=0.1;
    return CVect3((v1.j-v0.j)/f, v0.i-v1.i, 0)/norm(v0);
}

//***************************  class CQuat  *********************************/
CQuat::CQuat(void)
{
}

CQuat::CQuat(double qq0, const CVect3 &qqv)
{
	q0=qq0, q1=qqv.i, q2=qqv.j, q3=qqv.k;
}

CQuat::CQuat(double qq0, double qq1, double qq2, double qq3)
{
	q0=qq0, q1=qq1, q2=qq2, q3=qq3;
}

CQuat::CQuat(const double *pdata)
{
	q0=*pdata++, q1=*pdata++, q2=*pdata++, q3=*pdata++;
}

CQuat CQuat::operator+(const CVect3 &phi) const
{
	CQuat qtmp = rv2q(-phi);
	return qtmp*(*this);
}

CQuat CQuat::operator-(const CVect3 &phi) const
{
	CQuat qtmp = rv2q(phi);
	return qtmp*(*this);
}

CVect3 CQuat::operator-(const CQuat &quat) const
{
	CQuat dq;
	
	dq = quat*(~(*this));
	if(dq.q0<0)
	{
		dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3;
	}
	double n2 = acos(dq.q0), f;
	if( sign(n2)!=0 )
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

CQuat CQuat::operator*(const CQuat &quat) const
{
	CQuat qtmp;
	qtmp.q0 = q0*quat.q0 - q1*quat.q1 - q2*quat.q2 - q3*quat.q3;
	qtmp.q1 = q0*quat.q1 + q1*quat.q0 + q2*quat.q3 - q3*quat.q2;
	qtmp.q2 = q0*quat.q2 + q2*quat.q0 + q3*quat.q1 - q1*quat.q3;
	qtmp.q3 = q0*quat.q3 + q3*quat.q0 + q1*quat.q2 - q2*quat.q1;
	return qtmp;
}

CQuat& CQuat::operator*=(const CQuat &quat)
{
	return (*this=*this*quat);
}

CQuat& CQuat::operator-=(const CVect3 &phi)
{
	CQuat qtmp = rv2q(phi);
	return (*this=qtmp*(*this));
}

CQuat operator~(const CQuat &q)
{
	return CQuat(q.q0,-q.q1,-q.q2,-q.q3);
}

CVect3 CQuat::operator*(const CVect3 &v) const
{
	CQuat qtmp;
	CVect3 vtmp;
	qtmp.q0 =         - q1*v.i - q2*v.j - q3*v.k;
	qtmp.q1 = q0*v.i           + q2*v.k - q3*v.j;
	qtmp.q2 = q0*v.j           + q3*v.i - q1*v.k;
	qtmp.q3 = q0*v.k           + q1*v.j - q2*v.i;
	vtmp.i = -qtmp.q0*q1 + qtmp.q1*q0 - qtmp.q2*q3 + qtmp.q3*q2;
	vtmp.j = -qtmp.q0*q2 + qtmp.q2*q0 - qtmp.q3*q1 + qtmp.q1*q3;
	vtmp.k = -qtmp.q0*q3 + qtmp.q3*q0 - qtmp.q1*q2 + qtmp.q2*q1;
	return vtmp;
}

void CQuat::SetYaw(double yaw)
{
	CVect3 att = q2att(*this);
	att.k = yaw;
	*this = a2qua(att);
}

double normlize(CQuat *q)
{
	double nq=sqrt(q->q0*q->q0+q->q1*q->q1+q->q2*q->q2+q->q3*q->q3);
	q->q0 /= nq, q->q1 /= nq, q->q2 /= nq, q->q3 /= nq;
	return nq;
}

CVect3 q2rv(const CQuat &q)
{
	CQuat dq;
	dq = q;
	if(dq.q0<0)  { dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3; }
	if(dq.q0>1.0) dq.q0=1.0;
	double n2 = acos(dq.q0), f;
	if(n2>1.0e-20)
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

CVect3 qq2phi(const CQuat &qcalcu, const CQuat &qreal)
{
    return q2rv(qreal*(~qcalcu));
}

CQuat addmu(const CQuat &q, const CVect3 &mu)
{
	return q*rv2q(mu);
}

CQuat UpDown(const CQuat &q)
{
	CVect3 att = q2att(q);
	att.i = -att.i; att.j += PI;
	return a2qua(att);
}

//***************************  class CMat3  *********************************/
CMat3::CMat3(void)
{
}

CMat3::CMat3(double xyz)
{
	e00=e01=e02 =e10=e11=e12 =e20=e21=e22 =xyz;
}

CMat3::CMat3(const double *pxyz)
{
	e00=*pxyz++,e01=*pxyz++,e02=*pxyz++,
	e10=*pxyz++,e11=*pxyz++,e12=*pxyz++,
	e20=*pxyz++,e21=*pxyz++,e22=*pxyz  ;
}

CMat3::CMat3(const float *pxyz)
{
	e00=*pxyz++,e01=*pxyz++,e02=*pxyz++,
	e10=*pxyz++,e11=*pxyz++,e12=*pxyz++,
	e20=*pxyz++,e21=*pxyz++,e22=*pxyz  ;
}

CMat3::CMat3(double xx, double yy, double zz)
{
	e00=xx, e11=yy, e22=zz;
	e01=e02 =e10=e12 =e20=e21 =0.0;
}

CMat3::CMat3(double xx, double xy, double xz, 
		  double yx, double yy, double yz,
		  double zx, double zy, double zz )
{
	e00=xx,e01=xy,e02=xz; e10=yx,e11=yy,e12=yz; e20=zx,e21=zy,e22=zz;
}

CMat3::CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2, BOOL isrow)
{
	if(isrow) {
		e00 = v0.i, e01 = v0.j, e02 = v0.k;
		e10 = v1.i, e11 = v1.j, e12 = v1.k;
		e20 = v2.i, e21 = v2.j, e22 = v2.k;
	}
	else {
		e00 = v0.i, e01 = v1.i, e02 = v2.i;
		e10 = v0.j, e11 = v1.j, e12 = v2.j;
		e20 = v0.k, e21 = v1.k, e22 = v2.k;
	}
}

CVect3 sv2att(const CVect3 &fb, double yaw0, const CVect3 &fn)
{
    CVect3 phi = fb*fn;
    double afa = acos(dot(fn,fb)/norm(fn)/norm(fb)), nphi = norm(phi);
	CVect3 att = q2att(rv2q(phi*(nphi<1.0e-10?1.0:afa/nphi)));   att.k = yaw0;
	return att;  // return C^n_b
}

CVect3 dv2att(const CVect3 &va1, const CVect3 &va2, const CVect3 &vb1, const CVect3 &vb2)
{
	CVect3 a=va1*va2, b=vb1*vb2, aa=a*va1, bb=b*vb1;
	if(IsZero(va1)||IsZero(a)||IsZero(aa)||IsZero(vb1)||IsZero(b)||IsZero(bb)) return O31;
	CMat3 Ma(va1/norm(va1),a/norm(a),aa/norm(aa)), Mb(vb1/norm(vb1),b/norm(b),bb/norm(bb));
	return m2att((~Ma)*(Mb));  // return C^a_b -> att
}

CVect3 mv2att(int n, const CVect3 *vai, ...)
{
	psinsassert(n>=2);
	CMat3 A(0.0);
	CVect3 *vbi;
	va_list vl;
	va_start(vl, vai);  vbi = va_arg(vl, CVect3*);  
	for(int i=0; i<n; i++)
	{ A += vxv(*vai,*vbi);  vai = va_arg(vl, CVect3*);  vbi = va_arg(vl, CVect3*); }
	va_end(vl);
	return m2att(sfoam(A));  // return C^a_b -> att
}

CVect3 vn2att(const CVect3 &vn)
{
	double vel = normXY(vn);
	if(vel<1.0e-6) return O31;
    return CVect3(atan2(vn.k, vel), 0, atan2(-vn.i, vn.j));
}

CVect3 atss(CVect3 &att, CVect3 &vn)
{
	CVect3 att1 = vn2att(vn);
	CVect3 as(att1.i-att.i, 0.0, att1.k-att.k);
	if(as.k>PI) as.k-=_2PI; else if(as.k<-PI) as.k+=_2PI;
	return as;
}

CMat3 operator-(const CMat3 &m)
{
	return CMat3(-m.e00,-m.e01,-m.e02,-m.e10,-m.e11,-m.e12,-m.e20,-m.e21,-m.e22);
}

CMat3 operator~(const CMat3 &m)
{
	return CMat3(m.e00,m.e10,m.e20, m.e01,m.e11,m.e21, m.e02,m.e12,m.e22);
}

CMat3 CMat3::operator*(const CMat3 &mat) const
{
	CMat3 mtmp;
	mtmp.e00 = e00*mat.e00 + e01*mat.e10 + e02*mat.e20;
	mtmp.e01 = e00*mat.e01 + e01*mat.e11 + e02*mat.e21;
	mtmp.e02 = e00*mat.e02 + e01*mat.e12 + e02*mat.e22;
	mtmp.e10 = e10*mat.e00 + e11*mat.e10 + e12*mat.e20;
	mtmp.e11 = e10*mat.e01 + e11*mat.e11 + e12*mat.e21;
	mtmp.e12 = e10*mat.e02 + e11*mat.e12 + e12*mat.e22;
	mtmp.e20 = e20*mat.e00 + e21*mat.e10 + e22*mat.e20;
	mtmp.e21 = e20*mat.e01 + e21*mat.e11 + e22*mat.e21;
	mtmp.e22 = e20*mat.e02 + e21*mat.e12 + e22*mat.e22;
	return mtmp;
}

CMat3 CMat3::operator+(const CMat3 &mat) const
{
	CMat3 mtmp;
	mtmp.e00 = e00 + mat.e00;  mtmp.e01 = e01 + mat.e01;  mtmp.e02 = e02 + mat.e02;  
	mtmp.e10 = e10 + mat.e10;  mtmp.e11 = e11 + mat.e11;  mtmp.e12 = e12 + mat.e12;  
	mtmp.e20 = e20 + mat.e20;  mtmp.e21 = e21 + mat.e21;  mtmp.e22 = e22 + mat.e22;  
	return mtmp;
}

CMat3& CMat3::operator+=(const CMat3 &mat)
{
	this->e00 += mat.e00;  this->e01 += mat.e01;  this->e02 += mat.e02;  
	this->e10 += mat.e10;  this->e11 += mat.e11;  this->e12 += mat.e12;  
	this->e20 += mat.e20;  this->e21 += mat.e21;  this->e22 += mat.e22;  
	return *this;
}

CMat3 CMat3::operator+(const CVect3 &v) const
{
	CMat3 mtmp=*this;
	mtmp.e00 += v.i;  mtmp.e11 += v.j;  mtmp.e22 += v.k;
	return mtmp;
}

CMat3& CMat3::operator+=(const CVect3 &v)
{
	this->e00 += v.i;  this->e11 += v.j;  this->e22 += v.k;  
	return *this;
}

CMat3 CMat3::operator-(const CMat3 &mat) const
{
	CMat3 mtmp;
	mtmp.e00 = e00 - mat.e00;  mtmp.e01 = e01 - mat.e01;  mtmp.e02 = e02 - mat.e02;  
	mtmp.e10 = e10 - mat.e10;  mtmp.e11 = e11 - mat.e11;  mtmp.e12 = e12 - mat.e12;  
	mtmp.e20 = e20 - mat.e20;  mtmp.e21 = e21 - mat.e21;  mtmp.e22 = e22 - mat.e22;  
	return mtmp;
}

CMat3 CMat3::operator*(double f) const
{
	return CMat3(e00*f,e01*f,e02*f, e10*f,e11*f,e12*f, e20*f,e21*f,e22*f);
}

double& CMat3::operator()(int r, int c)
{
	return (&this->e00)[r*3+c];
}

void CMat3::SetRow(int i, const CVect3 &v)
{
	double *p=&e00+i*3;
	*p=v.i, *(p+1)=v.j, *(p+2) = v.k;
}

void CMat3::SetClm(int i, const CVect3 &v)
{
	double *p=&e00+i;
	*p=v.i, *(p+3)=v.j, *(p+6) = v.k;
}

CVect3 CMat3::GetRow(int i) const
{
	const double *p=&e00+i*3;
	return CVect3(*p,*(p+1),*(p+2));
}

CVect3 CMat3::GetClm(int i) const
{
	const double *p=&e00+i;
	return CVect3(*p,*(p+3),*(p+6));
}

CMat3 Rot(double angle, char axis)
{
	double s=sin(angle), c=cos(angle);
	if(axis=='x'||axis=='X')		return CMat3(1,  0, 0,   0, c, -s,   0, s, c);
	else if(axis=='y'||axis=='Y')	return CMat3(c,  0, s,   0, 1,  0,  -s, 0, c);
	else							return CMat3(c, -s, 0,   s, c,  0,   0, 0, 1);
}

CVect3 rotz(const CVect3 &v, double angle)
{
	double s=sin(angle), c=cos(angle);
	return CVect3(v.i*c-v.j*s, v.i*s+v.j*c, v.k);
}

CMat3 rcijk(const CMat3 &m, int ijk)
{
	switch(ijk)
	{
//	case 012: return m; break;
	case 021: return CMat3(m.e00,m.e02,m.e01, m.e20,m.e22,m.e21, m.e10,m.e12,m.e11); break;
	case 102: return CMat3(m.e11,m.e10,m.e12, m.e01,m.e00,m.e02, m.e21,m.e20,m.e22); break;
	case 120: return CMat3(m.e11,m.e12,m.e10, m.e21,m.e22,m.e20, m.e01,m.e02,m.e00); break;
	case 201: return CMat3(m.e22,m.e20,m.e21, m.e02,m.e00,m.e01, m.e12,m.e10,m.e11); break;
	case 210: return CMat3(m.e22,m.e21,m.e20, m.e12,m.e11,m.e10, m.e02,m.e01,m.e00); break;
	}
	return m;
}

void symmetry(CMat3 &m)
{
	m.e01 = m.e10 = (m.e01+m.e10)*0.5;
	m.e02 = m.e20 = (m.e02+m.e20)*0.5;
	m.e12 = m.e21 = (m.e12+m.e21)*0.5;
}

CMat3 operator*(double f, const CMat3 &m)
{
	return CMat3(m.e00*f,m.e01*f,m.e02*f, m.e10*f,m.e11*f,m.e12*f, m.e20*f,m.e21*f,m.e22*f);
}

CVect3 CMat3::operator*(const CVect3 &v) const
{
	return CVect3(e00*v.i+e01*v.j+e02*v.k,e10*v.i+e11*v.j+e12*v.k,e20*v.i+e21*v.j+e22*v.k);
	//stacksize();
}

CMat3 dotmul(const CMat3 &m1, const CMat3 &m2)
{
	CMat3 m;
	m.e00 = m1.e00*m2.e00, m.e01 = m1.e01*m2.e01, m.e02 = m1.e02*m2.e02; 
	m.e10 = m1.e10*m2.e10, m.e11 = m1.e11*m2.e11, m.e12 = m1.e12*m2.e12; 
	m.e20 = m1.e20*m2.e20, m.e21 = m1.e21*m2.e21, m.e22 = m1.e22*m2.e22; 
	return m;
}

double det(const CMat3 &m)
{
	return m.e00*(m.e11*m.e22-m.e12*m.e21) - m.e01*(m.e10*m.e22-m.e12*m.e20) + m.e02*(m.e10*m.e21-m.e11*m.e20);
}

double trace(const CMat3 &m)
{
	return (m.e00+m.e11+m.e22);
}

CMat3 pow(const CMat3 &m, int k)
{
	CMat3 mm = m;
	for(int i=1; i<k; i++)	mm = mm*m;
	return mm;
}

CQuat a2qua(double pitch, double roll, double yaw)
{
	pitch /= 2.0, roll /= 2.0, yaw /= 2.0;
    double	sp = sin(pitch), sr = sin(roll), sy = sin(yaw), 
			cp = cos(pitch), cr = cos(roll), cy = cos(yaw);
	CQuat qnb;
    qnb.q0 = cp*cr*cy - sp*sr*sy;
    qnb.q1 = sp*cr*cy - cp*sr*sy;
    qnb.q2 = cp*sr*cy + sp*cr*sy;
    qnb.q3 = cp*cr*sy + sp*sr*cy;
	return qnb;
}

CQuat a2qua(const CVect3 &att)
{
	return a2qua(att.i, att.j, att.k);
}

CMat3 a2mat(const CVect3 &att)
{
	double	si = sin(att.i), ci = cos(att.i),
			sj = sin(att.j), cj = cos(att.j),
			sk = sin(att.k), ck = cos(att.k);
	CMat3 Cnb;
	Cnb.e00 =  cj*ck - si*sj*sk;	Cnb.e01 =  -ci*sk;	Cnb.e02 = sj*ck + si*cj*sk;
	Cnb.e10 =  cj*sk + si*sj*ck;	Cnb.e11 =  ci*ck;	Cnb.e12 = sj*sk - si*cj*ck;
	Cnb.e20 = -ci*sj;				Cnb.e21 =  si;		Cnb.e22 = ci*cj;
	return Cnb;
}

CMat3 ar2mat(const CVect3 &attr)  // reversed Euler angles to DCM
{
	double	si = sin(attr.i), ci = cos(attr.i),
			sj = sin(attr.j), cj = cos(attr.j),
			sk = sin(attr.k), ck = cos(attr.k);
	CMat3 Cnb;
	Cnb.e00 =  cj*ck;	Cnb.e01 =  si*sj*ck-ci*sk;	Cnb.e02 = ci*sj*ck + si*sk;
	Cnb.e10 =  cj*sk;	Cnb.e11 =  si*sj*sk+ci*ck;	Cnb.e12 = ci*sj*sk - si*ck;
	Cnb.e20 = -sj;		Cnb.e21 =  si*cj;			Cnb.e22 = ci*cj;
	return Cnb;
}

CQuat ar2qua(const CVect3 &attr)
{
	return m2qua(ar2mat(attr));
}

CVect3 m2att(const CMat3 &Cnb)
{
	CVect3 att;
	att.i = asinEx(Cnb.e21);
	att.j = atan2Ex(-Cnb.e20, Cnb.e22);
	att.k = atan2Ex(-Cnb.e01, Cnb.e11);
	return att;
}

CVect3 m2attr(const CMat3 &Cnb)
{
	CVect3 attr;
	attr.i = atan2Ex(Cnb.e21, Cnb.e22);
	attr.j = asinEx(-Cnb.e20);
	attr.k = atan2Ex(Cnb.e10, Cnb.e00);
	return attr;
}

CVect3 q2attr(const CQuat &qnb)
{
	return m2attr(q2mat(qnb));
}

CVect3 m2rv(const CMat3 &Cnb)
{
	double phi=acos((trace(Cnb)-1.0)/2), afa;
	afa = (-1.0e-10<phi&&phi<1.0e-10) ? 1.0/2 : phi/(2*sin(phi));
	return CVect3(Cnb.e21-Cnb.e12, Cnb.e02-Cnb.e20, Cnb.e10-Cnb.e01)*afa;
}

CQuat m2qua(const CMat3 &Cnb)
{
	double q0, q1, q2, q3, qq4;
    if(Cnb.e00>=Cnb.e11+Cnb.e22)
	{
        q1 = 0.5*sqrt(1+Cnb.e00-Cnb.e11-Cnb.e22);  qq4 = 4*q1;
        q0 = (Cnb.e21-Cnb.e12)/qq4; q2 = (Cnb.e01+Cnb.e10)/qq4; q3 = (Cnb.e02+Cnb.e20)/qq4;
	}
    else if(Cnb.e11>=Cnb.e00+Cnb.e22)
	{
        q2 = 0.5*sqrt(1-Cnb.e00+Cnb.e11-Cnb.e22);  qq4 = 4*q2;
        q0 = (Cnb.e02-Cnb.e20)/qq4; q1 = (Cnb.e01+Cnb.e10)/qq4; q3 = (Cnb.e12+Cnb.e21)/qq4;
	}
    else if(Cnb.e22>=Cnb.e00+Cnb.e11)
	{
        q3 = 0.5*sqrt(1-Cnb.e00-Cnb.e11+Cnb.e22);  qq4 = 4*q3;
        q0 = (Cnb.e10-Cnb.e01)/qq4; q1 = (Cnb.e02+Cnb.e20)/qq4; q2 = (Cnb.e12+Cnb.e21)/qq4;
	}
    else
	{
        q0 = 0.5*sqrt(1+Cnb.e00+Cnb.e11+Cnb.e22);  qq4 = 4*q0;
        q1 = (Cnb.e21-Cnb.e12)/qq4; q2 = (Cnb.e02-Cnb.e20)/qq4; q3 = (Cnb.e10-Cnb.e01)/qq4;
	}
	double nq = sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
	q0 /= nq; q1 /= nq; q2 /= nq; q3 /= nq;
	return CQuat(q0, q1, q2, q3);
}

CVect3 q2att(const CQuat &qnb)
{
	double	q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
			q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
			q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
			q44 = qnb.q3*qnb.q3;
	CVect3 att;
	att.i = asinEx(2*(q34+q12));
	att.j = atan2Ex(-2*(q24-q13), q11-q22-q33+q44);
	att.k = atan2Ex(-2*(q23-q14), q11-q22+q33-q44);
	return att;
}

CMat3 q2mat(const CQuat &qnb)
{
	double	q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
			q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
			q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
			q44 = qnb.q3*qnb.q3;
	CMat3 Cnb;
    Cnb.e00 = q11+q22-q33-q44,  Cnb.e01 = 2*(q23-q14),     Cnb.e02 = 2*(q24+q13),
	Cnb.e10 = 2*(q23+q14),      Cnb.e11 = q11-q22+q33-q44, Cnb.e12 = 2*(q34-q12),
	Cnb.e20 = 2*(q24-q13),      Cnb.e21 = 2*(q34+q12),     Cnb.e22 = q11-q22-q33+q44;
	return Cnb;
}

CMat3 Ka2Cba(const CMat3 &Ka, CVect3 &Sfa)
{
    CMat3 iKa = inv(Ka);
    Sfa = CVect3(norm(*(CVect3*)&iKa.e00), norm(*(CVect3*)&iKa.e10), norm(*(CVect3*)&iKa.e20));
    CMat3 Cba = ~(inv(diag(Sfa))*iKa);
	return Cba;
}

CMat3 Cba2Ka(const CMat3 &Cba, const CVect3 &Sfa)
{
   CMat3 Ka = inv(diag(Sfa)*(~Cba));
   return Ka;
}

void Ka22Kpn(const CVect3 &Ka1, const CVect3 &Ka2, CVect3 &Kap, CVect3 &Kan)
{
	CVect3 dKa1=G0*Ka2;
	Kap = dotmul(Ka1,One31+dKa1), Kan = dotmul(Ka1,One31-dKa1);
}

void Kpn2Ka2(const CVect3 &Kap, const CVect3 &Kan, CVect3 &Ka1, CVect3 &Ka2)
{
	Ka1 = (Kap+Kan)/2.0;
	Ka2 = (dotdiv(Kap,Ka1)-One31)*(1/G0);
}

void KgMdf(CMat3 &Kg, const double *dKg, int typ)
{
//	CMat3 IdKG = I33-(~(*(CMat3*)&FBXk.dd[12]));
//	sins.imu.Kg = IdKG*sins.imu.Kg;			// dKG
	CMat3 Cbg0, Cbg; CVect3 Sfg0, Sfg;
	if(typ) Cbg0 = Ka2Cba(Kg, Sfg0);
	double de00=1.0-dKg[0], de11=1.0-dKg[4], de22=1.0-dKg[8], e00, e01, e02, e10, e11, e12;
	   e00 =  de00  *Kg.e00 - dKg[3]*Kg.e10 - dKg[6]*Kg.e20;
	   e01 =  de00  *Kg.e01 - dKg[3]*Kg.e11 - dKg[6]*Kg.e21;
	   e02 =  de00  *Kg.e02 - dKg[3]*Kg.e12 - dKg[6]*Kg.e22;
	   e10 = -dKg[1]*Kg.e00 + de11  *Kg.e10 - dKg[7]*Kg.e20;
	   e11 = -dKg[1]*Kg.e01 + de11  *Kg.e11 - dKg[7]*Kg.e21;
	   e12 = -dKg[1]*Kg.e02 + de11  *Kg.e12 - dKg[7]*Kg.e22;
	Kg.e20 = -dKg[2]*Kg.e00 - dKg[5]*Kg.e10 + de22  *Kg.e20;
	Kg.e21 = -dKg[2]*Kg.e01 - dKg[5]*Kg.e11 + de22  *Kg.e21;
	Kg.e22 = -dKg[2]*Kg.e02 - dKg[5]*Kg.e12 + de22  *Kg.e22;
	Kg.e00=e00, Kg.e01=e01, Kg.e02=e02, Kg.e10=e10, Kg.e11=e11, Kg.e12=e12;
	if(typ) {  // bad
		Cbg = Ka2Cba(Kg, Sfg);
		Kg = typ==1 ? Cba2Ka(Cbg, Sfg0) : Cba2Ka(Cbg0, Sfg);
	}
}

void KaMdf(CMat3 &Ka, const double *dKa, int typ)
{
//	CMat3 IdKA(1.0-FBXk.dd[21],            0.0,           0.0,
//			      -FBXk.dd[22],1.0-FBXk.dd[24],           0.0,
//				  -FBXk.dd[23],   -FBXk.dd[25],1.0-FBXk.dd[26]);
//	sins.imu.Ka = IdKA*sins.imu.Ka;			// dKA
	CMat3 Cba0, Cba; CVect3 Sfa0, Sfa;
	if(typ) Cba0 = Ka2Cba(Ka, Sfa0);
	double de00=1.0-dKa[0], de11=1.0-dKa[3], de22=1.0-dKa[5];
	Ka.e20 = -dKa[2]*Ka.e00 - dKa[4]*Ka.e10 + de22  *Ka.e20;
	Ka.e21 = -dKa[2]*Ka.e01 - dKa[4]*Ka.e11 + de22  *Ka.e21;
	Ka.e22 = -dKa[2]*Ka.e02 - dKa[4]*Ka.e12 + de22  *Ka.e22;
	Ka.e10 = -dKa[1]*Ka.e00 + de11  *Ka.e10;
	Ka.e11 = -dKa[1]*Ka.e01 + de11  *Ka.e11;
	Ka.e12 = -dKa[1]*Ka.e02 + de11  *Ka.e12;
	Ka.e00 =  de00  *Ka.e00;
	Ka.e01 =  de00  *Ka.e01;
	Ka.e02 =  de00  *Ka.e02;
	if(typ) {
		Cba = Ka2Cba(Ka, Sfa);
		Ka = typ==1 ? Cba2Ka(Cba, Sfa0) : Cba2Ka(Cba0, Sfa);
	}
}

CMat3 adj(const CMat3 &m)
{
	CMat3 mtmp;
	mtmp.e00 =  (m.e11*m.e22-m.e12*m.e21);
	mtmp.e10 = -(m.e10*m.e22-m.e12*m.e20);
	mtmp.e20 =  (m.e10*m.e21-m.e11*m.e20);
	mtmp.e01 = -(m.e01*m.e22-m.e02*m.e21);
	mtmp.e11 =  (m.e00*m.e22-m.e02*m.e20);
	mtmp.e21 = -(m.e00*m.e21-m.e01*m.e20);
	mtmp.e02 =  (m.e01*m.e12-m.e02*m.e11);
	mtmp.e12 = -(m.e00*m.e12-m.e02*m.e10);
	mtmp.e22 =  (m.e00*m.e11-m.e01*m.e10);
	return mtmp;
}

CMat3 inv(const CMat3 &m)
{
	CMat3 adjm = adj(m);
	double detm = m.e00*adjm.e00 + m.e01*adjm.e10 + m.e02*adjm.e20;
	return adjm*(1.0/detm);
}

CVect3 diag(const CMat3 &m)
{
	return CVect3(m.e00, m.e11, m.e22);
}

CMat3 diag(const CVect3 &v)
{
	return CMat3(v.i,0,0, 0,v.j,0, 0,0,v.k);
}

CMat3 diag(double ii, double jj, double kk)
{
	return CMat3(ii,0,0, 0,jj,0, 0,0,kk);
}

CMat3 askew(const CMat3 &m, int I)
{
	CMat3 m1;
	m1.e01=(m.e01-m.e10)/2;  m1.e02=(m.e02-m.e20)/2;  m1.e12=(m.e12-m.e21)/2;
	m1.e10=-m1.e01; m1.e20=-m1.e02; m1.e21=-m1.e12;
	if(I==0)  m1.e00=m1.e11=m1.e22=0.0;
	else if(I==1)  m1.e00=m1.e11=m1.e22=1.0;
	else  { m1.e00=m.e00, m1.e11=m.e11, m1.e22=m.e22; }
	return m1;
}

CMat3 MMT(const CMat3 &m1, const CMat3 &m2)
{
	CMat3 mtmp; const CMat3 *pm2;
	pm2 = (&m2==&I33) ? &m1 : &m2;
	mtmp.e00 = m1.e00*pm2->e00 + m1.e01*pm2->e01 + m1.e02*pm2->e02;
	mtmp.e01 = m1.e00*pm2->e10 + m1.e01*pm2->e11 + m1.e02*pm2->e12;
	mtmp.e02 = m1.e00*pm2->e20 + m1.e01*pm2->e21 + m1.e02*pm2->e22;
	mtmp.e11 = m1.e10*pm2->e10 + m1.e11*pm2->e11 + m1.e12*pm2->e12;
	mtmp.e12 = m1.e10*pm2->e20 + m1.e11*pm2->e21 + m1.e12*pm2->e22;
	mtmp.e22 = m1.e20*pm2->e20 + m1.e21*pm2->e21 + m1.e22*pm2->e22;
	mtmp.e10 = mtmp.e01;
	mtmp.e20 = mtmp.e02;
	mtmp.e21 = mtmp.e12;
	return mtmp;
}

double trMMT(const CMat3 &m1, const CMat3 &m2)
{
	const CMat3 *pm2;
	pm2 = (&m2==&I33) ? &m1 : &m2;
	return	m1.e00*pm2->e00 + m1.e01*pm2->e01 + m1.e02*pm2->e02 +
			m1.e10*pm2->e10 + m1.e11*pm2->e11 + m1.e12*pm2->e12 +
			m1.e20*pm2->e20 + m1.e21*pm2->e21 + m1.e22*pm2->e22;
}

double norm(const CMat3 &m)
{
	return sqrt(trMMT(m));
}

static double maxrt4(double b, double c, double d) // max real root for x^4+bx^2+cx+d=0
{
	double D=-8*b, E=-8*c, F=16*b*b-64*d, E2=E*E,
		A=D*D-3*F, B=D*F-9*E2;  //C=F*F-3*D*E2;// Delta=B*B-4*A*C;
	double sA, y, theta3, sAct3, sAst3, y1, y2, y3, sy1, sy2, sy3;
	if(A<0.0) { sA = 0.0,      y = 0.0; }
	else      { sA = sqrt(A),  y = (1.5*B/A-D)/sA; }  // y=(3*B-2*A*D)/(2*A*sA);
	if(y>1.0) y=1.0; else if(y<-1.0) y=-1.0;
	theta3 = acos(y)/3;
	sAct3 = sA*cos(theta3); sAst3 = sA*sqrt3*sin(theta3);
	y1 = D-2*sAct3      ;  sy1 = y1<0.0 ? 0.0 : sqrt(y1);
	y2 = D+  sAct3+sAst3;  sy2 = y2<0.0 ? 0.0 : sqrt(y2);
	y3 = D+  sAct3-sAst3;  sy3 = y3<0.0 ? 0.0 : sqrt(y3);
	if(E<0.0) sy1=-sy1;
	return (sy1+sy2+sy3)/(4*sqrt3);
}

CMat3 sfoam(const CMat3 &B, int iter)
{
    CMat3 adjBp=adj(~B), BBp=MMT(B), C;
    double detB=det(B), adjBp2=trMMT(adjBp), B2=trace(BBp), cnd=0.0,
		lambda=0.0, lambda2=0.0, kappa=0.0, zeta=0.0, Psi=0.0, dPsi=0.0, dlambda=0.0;
	if(B2<EPS||adjBp2<EPS) return I33;  // rank(B)<2
	cnd = detB*detB/B2/adjBp2;  // 1/cond(B)^2
	lambda = cnd>1.0e-12 ? maxrt4(-2*B2,-8*detB,B2*B2-4*adjBp2) : sqrt(B2+2*sqrt(adjBp2));
    for(int k=1; k<iter; k++)
	{
		lambda2 = lambda*lambda;
        kappa = (lambda2-B2)/2.0;
        zeta = kappa*lambda - detB;
		if(zeta<EPS) return I33;  // singular value s2+s3=0
		Psi = (lambda2-B2); Psi = Psi*Psi - 8.0*lambda*detB - 4.0*adjBp2;
		dPsi = 8.0*zeta;
		dlambda = Psi / dPsi;
        lambda = lambda - dlambda;
		if(fabs(dlambda/lambda)<1.0e-15) break;
    }
    C = ((kappa+B2)*B+lambda*adjBp-BBp*B) * (1.0/zeta);
	double nC=trMMT(C);
    return (nC<3.0-1.0e-6||nC>3.0+1.0e-6) ? I33 : C;
}

//***************************  class CMat  *********************************/
CMat::CMat(void)
{
#ifdef PSINS_MAT_COUNT
	#pragma message("  PSINS_MAT_COUNT")
	if(iMax<++iCount) iMax = iCount;
#endif
}
	
CMat::CMat(int row0, int clm0)
{
#ifdef PSINS_MAT_COUNT
	if(iMax<++iCount) iMax = iCount;
#endif
	row=row0; clm=clm0; rc=row*clm;
}

CMat::CMat(int row0, int clm0, double f)
{
#ifdef PSINS_MAT_COUNT
	if(iMax<++iCount) iMax = iCount;
#endif
	row=row0; clm=clm0; rc=row*clm;
	for(double *pd=dd, *pEnd=&dd[rc]; pd<pEnd; pd++)  *pd = f;
}

CMat::CMat(int row0, int clm0, double f, double f1, ...)
{
#ifdef PSINS_MAT_COUNT
	if(iMax<++iCount) iMax = iCount;
#endif
	row=row0; clm=clm0; rc=row*clm;
	va_list vl;
	dd[0] = f;
	va_start(vl, f1);
	for(int i=1; i<rc; i++)
	{ if(f>2*INF) break;  f = va_arg(vl, double);  dd[i] = f;}
	va_end(vl);
}

CMat::CMat(int row0, int clm0, const double *pf)
{
#ifdef PSINS_MAT_COUNT
	if(iMax<++iCount) iMax = iCount;
#endif
	row=row0; clm=clm0; rc=row*clm;
	memcpy(dd, pf, rc*sizeof(double));
}

CMat::CMat(int clm0, const CVect *pv, ...)
{
//	this->CMat(pv->row, clm0);
	row=pv->row; clm=clm0; rc=row*clm;
	va_list vl;
	va_start(vl, pv);
	for(int i=0; i<clm; i++) {
		SetClm(i, *pv);
		pv = va_arg(vl, const CVect*);
	}
	va_end(vl);
}

#ifdef PSINS_MAT_COUNT
int CMat::iCount=0, CMat::iMax=0;
CMat::~CMat(void)
{
	iCount--;
}
#endif

void CMat::Reset(int row0, int clm0)
{
	row=row0; clm=clm0; rc=row*clm;
}

void CMat::Clear(void)
{
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)	*p = 0.0;
}

CMat CMat::operator*(const CMat &m0) const
{
#ifdef PSINS_MAT_COUNT
	++iCount;
#endif
	psinsassert(this->clm==m0.row);
	CMat mtmp(this->row,m0.clm);
	int m=this->row, k=this->clm, n=m0.clm;
	double *p=mtmp.dd; const double *p1i=this->dd, *p2=m0.dd;
	for(int i=0; i<m; i++,p1i+=k)
	{
		for(int j=0; j<n; j++)
		{
			double f=0.0; const double *p1is=p1i, *p1isEnd=&p1i[k], *p2sj=&p2[j];
			for(; p1is<p1isEnd; p1is++,p2sj+=n)
				f += (*p1is) * (*p2sj);
			*p++ = f;
		}
	}
	//stacksize();
	return mtmp;
}

CMat CMat::operator*(const CMat3 &m0) const
{
#ifdef PSINS_MAT_COUNT
	++iCount;
#endif
	psinsassert(this->clm==3);
	CMat mtmp(3,3);  mtmp.SetMat3(0,0,m0);
	return *this*mtmp;
}

CVect CMat::operator*(const CVect &v) const
{
	psinsassert(this->clm==v.row);
	CVect vtmp(this->row);
	double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row]; const double *p1ij=this->dd, *p2End=&v.dd[v.row];
	for(; p<pEnd; p++)
	{
		double f=0.0; const double *p2j=v.dd;
		for(; p2j<p2End; p1ij++,p2j++)	f += (*p1ij) * (*p2j);
		*p = f;
	}
	//stacksize();
	return vtmp;
}

CMat CMat::operator+(const CMat &m0) const
{
#ifdef PSINS_MAT_COUNT
	++iCount;
#endif
	psinsassert(row==m0.row&&clm==m0.clm);
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd, *p2=m0.dd;
	while(p<pEnd)
	{ *p++ = (*p1++) + (*p2++); } 
	return mtmp;
}

CMat& CMat::operator+=(const CVect &v)
{
	psinsassert(row==v.row||clm==v.clm);
	int row1 = row+1;
	double *p=dd, *pEnd=&dd[rc];
	for(const double *p1=v.dd; p<pEnd; p+=row1, p1++)	*p += *p1;
	return *this;
}

CMat CMat::operator-(const CMat &m0) const
{
#ifdef PSINS_MAT_COUNT
	++iCount;
#endif
	psinsassert(row==m0.row&&clm==m0.clm);
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd, *p2=m0.dd;
	while(p<pEnd)
	{ *p++ = (*p1++) - (*p2++); } 
	return mtmp;
}

CMat CMat::operator*(double f) const
{
#ifdef PSINS_MAT_COUNT
	++iCount;
#endif
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd;
	while(p<pEnd)
	{ *p++ = (*p1++) * f; } 
	return mtmp;
}

CMat& CMat::operator=(double f)
{
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)  { *p = f; }
	return *this;
}

CMat& CMat::operator+=(const CMat &m0)
{
	psinsassert(row==m0.row&&clm==m0.clm);
	double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
	while(p<pEnd)
	{ *p++ += *p1++; } 
	return *this;
}

CMat& CMat::operator-=(const CMat &m0)
{
	psinsassert(row==m0.row&&clm==m0.clm);
	double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
	while(p<pEnd)
	{ *p++ -= *p1++; } 
	//stacksize();
	return *this;
}

CMat& CMat::operator*=(double f)
{
	double *p=dd, *pEnd=&dd[rc];
	while(p<pEnd)
	{ *p++ *= f; } 
	return *this;
}

CMat& CMat::operator++()
{
	int row1=row+1;
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p+=row1)	*p += 1.0;
	return *this;
}

CMat operator~(const CMat &m0)
{
#ifdef PSINS_MAT_COUNT
	++CMat::iCount;
#endif
	CMat mtmp(m0.clm,m0.row);
	const double *pm=m0.dd;
	for(int i=0; i<m0.row; i++)
	{ for(int j=i; j<m0.rc; j+=m0.row) mtmp.dd[j] = *pm++; }
	return mtmp;
}

void symmetry(CMat &m)
{
	psinsassert(m.row==m.clm);
	double *prow0=&m.dd[1], *prowEnd=&m.dd[m.clm], *pclm0=&m.dd[m.clm], *pEnd=&m.dd[m.rc];
	for(int clm1=m.clm+1; prow0<pEnd; prow0+=clm1,pclm0+=clm1,prowEnd+=m.clm)
	{
		for(double *prow=prow0,*pclm=pclm0; prow<prowEnd; prow++,pclm+=m.clm)
			*prow = *pclm = (*prow+*pclm)*0.5;
	}
}

double trace(const CMat &m)
{
	psinsassert(m.row==m.clm);
	int row1 = m.row+1;
	double s = 0.0;
	for(const double *p=m.dd, *pEnd=&m.dd[m.rc]; p<pEnd; p+=row1)  s+=*p;
	return s;
}

CMat dotmul(const CMat &m1, const CMat &m2)
{
	psinsassert(m1.row==m2.row && m1.clm==m2.clm);
	CMat res(m1.row,m1.clm);
	double *p=res.dd;
	for(const double *p1=m1.dd, *p2=m2.dd, *pEnd=&m1.dd[m1.rc]; p1<pEnd; p1++,p2++,p++)  { *p = (*p1)*(*p2); }
	return res;
}

double& CMat::operator()(int r, int c)
{
	if(c<0) c = r;
	return this->dd[r*this->clm+c];
}

void CMat::SetRow(int i, double f, ...)
{
	va_list vl;
	va_start(vl, f);
	for(double *p=&dd[i*clm], *pEnd=p+clm; p<pEnd; p++)
	{ *p = f;  f = va_arg(vl, double);	}
	va_end(vl);
	return;
}

void CMat::SetRow(int i, const CVect &v)
{
	psinsassert(clm==v.clm);
	const double *p=v.dd;
	for(double *p1=&dd[i*clm],*pEnd=p1+clm; p1<pEnd; p++,p1++) *p1 = *p;
	return;
}

void CMat::SetClm(int j, double f, ...)
{
	va_list vl;
	va_start(vl, f);
	for(double *p=&dd[j], *pEnd=&p[rc]; p<pEnd; p+=clm)
	{ *p = f;  f = va_arg(vl, double);	}
	va_end(vl);
	return;
}

void CMat::SetClm(int j, const CVect &v)
{
	psinsassert(row==v.row);
	const double *p=v.dd;
	for(double *p1=&dd[j],*pEnd=&dd[rc]; p1<pEnd; p++,p1+=clm) *p1 = *p;
	return;
}

void CMat::SetClmVect3(int i, int j, const CVect3 &v)
{
	double *p=&dd[i*clm+j];
	*p = v.i; p += clm;
	*p = v.j; p += clm;
	*p = v.k;
}

void CMat::SetClmVect3(int i, int j, const CVect3 &v, const CVect3 &v1)
{
	double *p=&dd[i*clm+j];
	*p = v.i; *(p+1) = v1.i; p += clm;
	*p = v.j; *(p+1) = v1.j; p += clm;
	*p = v.k; *(p+1) = v1.k;
}

void CMat::SetClmVect3(int i, int j, const CVect3 &v, const CVect3 &v1, const CVect3 &v2)
{
	double *p=&dd[i*clm+j];
	*p = v.i; *(p+1) = v1.i; *(p+2) = v2.i; p += clm;
	*p = v.j; *(p+1) = v1.j; *(p+2) = v2.j; p += clm;
	*p = v.k; *(p+1) = v1.k; *(p+2) = v2.k; 
//	SetMat3(i,j, CMat3(v, v1, v2));
}

void CMat::SetRowVect3(int i, int j, const CVect3 &v)
{
	*(CVect3*)&dd[i*clm+j] = v;
}

void CMat::SetRowVect3(int i, int j, const CVect3 &v, const CVect3 &v1)
{
	CVect3 *p=(CVect3*)&dd[i*clm+j];
	*p++ = v;  *p = v1;
}

void CMat::SetRowVect3(int i, int j, const CVect3 &v, const CVect3 &v1, const CVect3 &v2)
{
	CVect3 *p=(CVect3*)&dd[i*clm+j];
	*p++ = v;  *p++ = v1;  *p = v2;
}

CVect3 CMat::GetRowVect3(int i, int j) const
{
	return *(CVect3*)&dd[i*clm+j];
}

CVect3 CMat::GetClmVect3(int i, int j) const
{
	CVect3 v;
	const double *p=&dd[i*clm+j];
	v.i = *p; p += clm;
	v.j = *p; p += clm;
	v.k = *p;
	return v;
}

void CMat::SetDiagVect3(int i, int j, const CVect3 &v)
{
	double *p=&dd[i*clm+j];
	*p = v.i;  p += clm+1;
	*p = v.j;  p += clm+1;
	*p = v.k;
}

CVect3 CMat::GetDiagVect3(int i, int j) const
{
	if(j==-1) j=i;
	CVect3 v;
	const double *p=&dd[i*clm+j];
	v.i = *p;  p += clm+1;
	v.j = *p;  p += clm+1;
	v.k = *p;
	return v;
}

void CMat::SetAskew(int i, int j, const CVect3 &v)
{
	double *p=&dd[i*clm+j];
	p[0] = 0.0; p[1] =-v.k; p[2] = v.j;  p += clm;
	p[0] = v.k; p[1] = 0.0; p[2] =-v.i;  p += clm;
	p[0] =-v.j; p[1] = v.i; p[2] = 0.0;
}

void CMat::SetMat3(int i, int j, const CMat3 &m)
{
	double *p=&dd[i*clm+j];
	*(CVect3*)p = *(CVect3*)&m.e00;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e10;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e20;
}

void CMat::SetMat3(int i, int j, const CMat3 &m, const CMat3 &m1)
{
	double *p=&dd[i*clm+j];
	*(CVect3*)p = *(CVect3*)&m.e00;  *(CVect3*)(p+3) = *(CVect3*)&m1.e00;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e10;  *(CVect3*)(p+3) = *(CVect3*)&m1.e10;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e20;  *(CVect3*)(p+3) = *(CVect3*)&m1.e20;  
}

void CMat::SetMat3(int i, int j, const CMat3 &m, const CMat3 &m1, const CMat3 &m2)
{
	double *p=&dd[i*clm+j];
	*(CVect3*)p = *(CVect3*)&m.e00;  *(CVect3*)(p+3) = *(CVect3*)&m1.e00;  *(CVect3*)(p+6) = *(CVect3*)&m2.e00;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e10;  *(CVect3*)(p+3) = *(CVect3*)&m1.e10;  *(CVect3*)(p+6) = *(CVect3*)&m2.e10;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e20;  *(CVect3*)(p+3) = *(CVect3*)&m1.e20;  *(CVect3*)(p+6) = *(CVect3*)&m2.e20;  
}

CMat3 CMat::GetMat3(int i, int j) const
{
	if(j==-1) j=i;
	CMat3 m;
	const double *p=&dd[i*clm+j];
	*(CVect3*)&m.e00 = *(CVect3*)p;  p += clm;
	*(CVect3*)&m.e10 = *(CVect3*)p;  p += clm;
	*(CVect3*)&m.e20 = *(CVect3*)p;
	return m;
}

void CMat::SubAddMat3(int i, int j, const CMat3 &m)
{
	double *p=&dd[i*clm+j];
	*(CVect3*)p += *(CVect3*)&m.e00;  p += clm;
	*(CVect3*)p += *(CVect3*)&m.e10;  p += clm;
	*(CVect3*)p += *(CVect3*)&m.e20;
}

CVect CMat::GetRow(int i) const
{
	CVect v(1, clm);
	const double *p1=&dd[i*clm], *pEnd=p1+clm;
	for(double *p=v.dd; p1<pEnd; p++,p1++) *p = *p1;
	return v;
}

void CMat::GetRow(CVect &v, int i)
{
	v.row=1, v.clm=v.rc=clm;
	const double *p1=&dd[i*clm], *pEnd=p1+clm;
	for(double *p=v.dd; p1<pEnd; p++,p1++) *p = *p1;
}

CVect CMat::GetClm(int j) const
{
	CVect v(row, 1);
	const double *p1=&dd[j], *pEnd=&dd[rc];
	for(double *p=v.dd; p1<pEnd; p++,p1+=clm) *p = *p1;
	return v;
}

CMat CMat::GetClm(int j1, ...) const
{
	va_list vl;
	va_start(vl, j1); int j=j1, n;
	for(n=0; n<MMD; n++) { 
		if(j<0) break;  
		j=va_arg(vl, int);
	}
	va_end(vl);
	//
	CMat m(row, n);
	va_start(vl, j1); j=j1;
	for(int i=0; i<n; i++) {
		const double *p1=&dd[j], *pEnd=&dd[rc];
		for(double *p=&m.dd[i]; p1<pEnd; p+=n,p1+=clm) *p = *p1;
		j=va_arg(vl, int);
	}
	va_end(vl);
	return m;
}

void CMat::GetClm(CVect &v, int j)
{
	v.row=v.rc=row, v.clm=1;
	const double *p1=&dd[j], *pEnd=&dd[rc];
	for(double *p=v.dd; p1<pEnd; p++,p1+=clm) *p = *p1;
}

void CMat::ZeroRow(int i)
{
	for(double *p=&dd[i*clm],*pEnd=p+clm; p<pEnd; p++) *p = 0.0;
	return;
}

void CMat::ZeroClm(int j)
{
	for(double *p=&dd[j],*pEnd=&dd[rc]; p<pEnd; p+=clm) *p = 0.0;
	return;
}

void CMat::SetDiag(double f, ...)
{
	*this = CMat(this->row, this->clm, 0.0);
	va_list vl;
	va_start(vl, f);
	double *p=dd, *pEnd=&dd[rc];
	for(int row1=row+1; p<pEnd; p+=row1)
	{ if(f>2*INF) break;  *p = f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CMat::SetDiag2(double f, ...)
{
	*this = CMat(this->row, this->clm, 0.0);
	va_list vl;
	va_start(vl, f);
	double *p=dd, *pEnd=&dd[rc];
	for(int row1=row+1; p<pEnd; p+=row1)
	{ if(f>2*INF) break;  *p = f*f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CMat::SetAscend(double f0, double df)
{
	for(int i=0; i<rc; i++, f0+=df)  dd[i]=f0;
}

double norm1(const CMat &m)
{
	return norm1(&m.dd[0], m.rc);
}

double normInf(const CMat &m)
{
	return normInf(&m.dd[0], m.rc);
}

CVect diag(const CMat &m)
{
	int row1 = m.row+1;
	CVect vtmp(m.row,1);
	double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row];
	for(const double *p1=m.dd; p<pEnd; p++, p1+=row1)	*p = *p1;
	return vtmp;
}

CMat eye(int n)
{
	CMat m(n,n, 0.0);
	double *p=m.dd, *pEnd=&m.dd[m.rc];
	for(n=n+1; p<pEnd; p+=n)  *p=1.0;
	return m;
}

CMat inv4(const CMat &m)
{
	psinsassert(m.clm==m.row&&m.clm==4);
	CMat3 A11(m.dd[0], m.dd[1], m.dd[2], m.dd[4], m.dd[5], m.dd[6], m.dd[8], m.dd[9], m.dd[10]), iA11;
	CVect3 A12(m.dd[3], m.dd[7], m.dd[11]), A21(m.dd[12], m.dd[13], m.dd[14]);
	double A22=m.dd[15], iA22;
	iA11 = inv(A11);  iA22 = 1.0/(A22-dot(A21*iA11,A12));
	CMat M(4,4);
	M.SetMat3(0,0, iA11+vxv(iA11*A12*iA22,A21)*iA11);  // by using matrix-inversion-lemma
	M.SetClmVect3(0,3, -iA11*A12*iA22);
	M.SetRowVect3(3,0, -iA22*A21*iA11);
	M.dd[15]=iA22;
	return M;
}

CMat inv6(const CMat &m)
{
	psinsassert(m.clm==m.row&&m.clm==6);
	CMat3 A11=m.GetMat3(0,0), iA11, A12=m.GetMat3(0,3), A21=m.GetMat3(3,0), A22=m.GetMat3(3,3), iA22;
	iA11 = inv(A11);  iA22 = inv(A22-A21*iA11*A12);
	CMat M(6,6);
	M.SetMat3(0,0, iA11+iA11*A12*iA22*A21*iA11);  // by using matrix-inversion-lemma
	M.SetMat3(0,3, -iA11*A12*iA22);
	M.SetMat3(3,0, -iA22*A21*iA11);
	M.SetMat3(3,3, iA22);
	return M;
}

CVect lss(const CMat &A, const CVect &y)
{
	CMat AT=~A;
	if(A.clm==4)
		return inv4(AT*A)*(AT*y);
	else if(A.clm==6)
		return inv6(AT*A)*(AT*y);
	else //error
		return Onen1;
}

void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r, int fast)
{
	psinsassert(m0.clm==m1.row);
	int rc0=r*m0.clm; fast=(r>=fast);
	double *p=&m.dd[rc0], *pEnd=p+m0.clm; const double *p0=&m0.dd[rc0], *p0End=p0+m0.clm, *p1j=m1.dd;
	for(; p<pEnd; p++,p1j++)
	{
		if(fast) { *p=p1j[r*m1.row]; continue; }
		double f=0.0; const double *p0j=p0, *p1jk=p1j;
		for(; p0j<p0End; p0j++,p1jk+=m1.clm)	 f += (*p0j) * (*p1jk);
		*p = f;
	}
}

void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r, int fast)
{
	psinsassert(m0.clm==m1.clm);
	int rc0=r*m0.clm, ifast=0;
	double *p=&m.dd[rc0], *pEnd=p+m0.clm; const double *p0=&m0.dd[rc0], *p0End=p0+m0.clm, *p1jk=m1.dd;
	for(; p<pEnd; p++,ifast++)
	{
		if(ifast>=fast) { *p=p0[ifast]; p1jk+=m1.clm; continue; }
		double f=0.0; const double *p0j=p0;
		for(; p0j<p0End; p0j++,p1jk++)	 f += (*p0j) * (*p1jk);
		*p = f;
	}
}

CMat diag(const CVect &v)
{
#ifdef PSINS_MAT_COUNT
	++CMat::iCount;
#endif
	int rc = v.row>v.clm ? v.row : v.clm, rc1=rc+1;
	CMat mtmp(rc,rc,0.0);
	double *p=mtmp.dd;
	for(const double *p1=v.dd, *p1End=&v.dd[rc]; p1<p1End; p+=rc1, p1++)	*p = *p1;
	return mtmp;
}

void DVMDVafa(const CVect &V, CMat &M, double afa)
{
	psinsassert(V.rc==M.row&&M.row==M.clm);
	int i = 0;
	const double *pv = V.dd;
	for(double vi=*pv, viafa=vi*afa; i<M.clm; i++,pv++,vi=*pv,viafa=vi*afa)
	{
		for(double *prow=&M.dd[i*M.clm],*prowEnd=prow+M.clm,*pclm=&M.dd[i]; prow<prowEnd; prow++,pclm+=M.row)
		{
			*prow *= vi;
			*pclm *= viafa;
		}
	}
}

//***************************  class CVect  *********************************/
CVect::CVect(void)
{
}

CVect::CVect(int row0, int clm0)
{
	if(clm0==1) { row=row0; clm=1;   }
	else		{ row=1;    clm=clm0;}
	rc = row*clm;
}

CVect::CVect(int row0, double f)
{
	row=row0; clm=1; rc=row*clm;
	row = row>PSINS_MATRIX_MAX_DIM?PSINS_MATRIX_MAX_DIM:row;
	for(int i=0;i<row;i++) dd[i]=f;
}

CVect::CVect(int row0, const double *pf)
{
	row=row0; clm=1; rc=row*clm;
	memcpy(dd, pf, row*sizeof(double));
}

CVect::CVect(int row0, double f, double f1, ...)
{
	row=row0; clm=1; rc=row*clm;
	psinsassert(row<=MMD&&clm<=MMD);
	dd[0] = f;
	va_list vl;
	va_start(vl, f1);
	rc=row>clm?row:clm;
	for(int i=1; i<rc; i++)
	{ if(f>2*INF) break;  dd[i] = f1;  f1 = va_arg(vl, double);	}
	va_end(vl);
}

CVect::CVect(const CVect3 &v)
{
	row=3; clm=1; rc=row*clm;
	dd[0]=v.i; dd[1]=v.j; dd[2]=v.k;
}

CVect::CVect(const CVect3 &v1, const CVect3 v2)
{
	row=6; clm=1; rc=row*clm;
	dd[0]=v1.i; dd[1]=v1.j; dd[2]=v1.k;
	dd[3]=v2.i; dd[4]=v2.j; dd[5]=v2.k;
}

void CVect::Clear(void)
{
	for(int i=0;i<rc;i++) dd[i]=0.0;
}

void CVect::Reset(int row0, int clm0)
{
	if(clm0==1) { row=row0; clm=1;   }
	else		{ row=1;    clm=clm0;}
	rc = row*clm;
}

CVect operator-(const CVect &v)
{
	CVect vtmp=v;
	for(double *p=&vtmp.dd[0], *pEnd=&vtmp.dd[vtmp.rc]; p<pEnd; p++) *p=-*p;
	return vtmp;
}

CVect operator~(const CVect &v)
{
	CVect vtmp=v;
	vtmp.row=v.clm; vtmp.clm=v.row;
	return vtmp;
}

CVect CVect::operator*(const CMat &m) const
{
	psinsassert(clm==m.row);
	CVect vtmp(row,clm);
	double *p=vtmp.dd; const double *p1End=&dd[clm];
	for(int j=0; j<clm; p++,j++)
	{
		double f=0.0; const double *p1j=dd, *p2jk=&m.dd[j];
		for(; p1j<p1End; p1j++,p2jk+=m.clm)	 f += (*p1j) * (*p2jk);
		*p = f;
	}
	return vtmp;
}

CMat CVect::operator*(const CVect &v) const
{
#ifdef MAT_STATISTIC
	++CMat::iCount;
#endif
	psinsassert(clm==v.row);
	CMat mtmp(row,v.clm);
	if(row==1 && v.clm==1)  // (1x1) = (1xn)*(nx1)
	{
		double f = 0.0;
		for(int i=0; i<clm; i++)  f += dd[i]*v.dd[i];
		mtmp.dd[0] = f;
	}
	else    // (nxn) = (nx1)*(1xn)
	{
		double *p=mtmp.dd;
		for(const double *p1=&dd[0],*p1End=&dd[rc],*p2End=&v.dd[rc]; p1<p1End; p1++)
		{
			for(const double *p2=&v.dd[0]; p2<p2End; p2++)  *p++ = *p1 * *p2;
		}
	}
	//stacksize();
	return mtmp;
}

CVect CVect::operator+(const CVect &v) const
{
	psinsassert(row==v.row&&clm==v.clm);
	const double *p2=v.dd, *p1=dd, *p1End=&dd[rc];
	CVect vtmp(row,clm);
	for(double *p=vtmp.dd; p1<p1End; p++,p1++,p2++)  { *p=*p1+*p2; }
	return vtmp;
}

CVect CVect::operator-(const CVect &v) const
{
	psinsassert(row==v.row&&clm==v.clm);
	const double *p2=v.dd, *p1=dd, *p1End=&dd[rc];
	CVect vtmp(row,clm);
	for(double *p=vtmp.dd; p1<p1End; p++,p1++,p2++)  { *p=*p1-*p2; }
	return vtmp;
}
	
CVect CVect::operator*(double f) const
{
	CVect vtmp(row,clm);
	const double *p1=dd,*p1End=&dd[rc];
	for(double *p=vtmp.dd; p1<p1End; p++,p1++)  { *p=*p1*f; }
	//stacksize();
	return vtmp;
}

CVect& CVect::operator=(double f)
{
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)  { *p = f; }
	return *this;
}

CVect& CVect::operator=(const double *pf)
{
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,pf++)  { *p = *pf; }
	return *this;
}

CVect& CVect::operator=(const CMat3 &m)
{
	row=9; clm=1; rc=9;
	memcpy(dd, &m.e00, 9*sizeof(double));
	return *this;
}

CVect& CVect::operator+=(const CVect &v)
{
	psinsassert(row==v.row&&clm==v.clm);
	const double *p1 = v.dd;
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,p1++)  { *p += *p1; }
	return *this;
}

CVect& CVect::operator-=(const CVect &v)
{
	psinsassert(row==v.row&&clm==v.clm);
	const double *p1 = v.dd;
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,p1++)  { *p -= *p1; }
	return *this;
}

CVect& CVect::operator*=(double f)
{
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)  { *p *= f; }
	return *this;
}

double dot(const CVect &v1, const CVect &v2)
{
	psinsassert(v1.row==v2.row && v1.clm==v2.clm);
	double res=0.0;
	for(const double *p1=v1.dd, *p2=v2.dd, *pEnd=&v1.dd[v1.rc]; p1<pEnd; p1++,p2++)  { res += (*p1)*(*p2); }
	return res;
}

CVect dotmul(const CVect &v1, const CVect &v2)
{
	psinsassert(v1.row==v2.row && v1.clm==v2.clm);
	CVect res(v1.row,v1.clm);
	double *p=res.dd;
	for(const double *p1=v1.dd, *p2=v2.dd, *pEnd=&v1.dd[v1.rc]; p1<pEnd; p1++,p2++,p++)  { *p = (*p1)*(*p2); }
	return res;
}

CVect pow(const CVect &v, int k)
{
	CVect pp = v;
	double *p, *pEnd=&pp.dd[pp.rc];
	for(int i=1; i<k; i++)
	{
		p=pp.dd;
		for(const double *p1=v.dd; p<pEnd; p++,p1++)
			*p *= *p1;
	}
	return pp;
}

CVect abs(const CVect &v)
{
	CVect res(v.row,v.clm);
	const double *p=v.dd, *pEnd=&v.dd[v.rc];
	for(double *p1=res.dd; p<pEnd; p++,p1++)  { *p1 = *p>0 ? *p : -*p; }
	return res;
}

double norm(const CVect &v)
{
	return norm(&v.dd[0], v.rc);
}

double norm1(const CVect &v)
{
	return norm1(&v.dd[0], v.rc);
}

double normInf(const CVect &v)
{
	return normInf(&v.dd[0], v.rc);
}

double& CVect::operator()(int r)
{
	return this->dd[r];
}

void CVect::Set(double f, ...)
{
	psinsassert(rc<=MMD);
	va_list vl;
	va_start(vl, f);
	for(int i=0; i<rc; i++)
	{ if(f>2*INF) break;  dd[i] = f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CVect::Set2(double f, ...)
{
	psinsassert(rc<=MMD);
	va_list vl;
	va_start(vl, f);
	for(int i=0; i<rc; i++)
	{ if(f>2*INF) break;  dd[i] = f*f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CVect::SetVect3(int i, const CVect3 &v)
{
	*(CVect3*)&dd[i] = v;
}

void CVect::Set2Vect3(int i, const CVect3 &v)
{
	dd[i++]=v.i*v.i; dd[i++]=v.j*v.j; dd[i]=v.k*v.k; 
}

void CVect::SetAscend(double f0, double df)
{
	for(int i=0; i<rc; i++, f0+=df)  dd[i]=f0;
}

void CVect::SetBit(unsigned int bit, double f)
{
	for(int i=0; i<rc; i++)		// assert(rc<32)
		if(bit&(0x01<<i)) dd[i]=f;
}

void CVect::SetBit(unsigned int bit, const CVect3 &v)
{
	const double *p=&v.i;
	for(int i=0; i<rc; i++)		// assert(rc<32)
		if(bit&(0x01<<i)) { dd[i]=*p++;  if(p>&v.k) p=&v.i; }
}

void CVect::Seti2j(int i, int j, double val)
{
	for(double *p=&dd[i], *pEnd=&dd[mmin(j,MMD-1)]; p<=pEnd; p++) *p = val;
}

CVect3 CVect::GetVect3(int i) const
{
	return *(CVect3*)&dd[i]; 
}

double mean(const CVect &v)
{
	double m=0.0;
	for(const double *p=v.dd, *pend=&v.dd[v.rc]; p<pend; p++) m += *p;
	return m/v.rc;
}

CVect sort(const CVect &v)
{
	CVect vtmp=v;
	double *pi=vtmp.dd, *pj, *pend=&vtmp.dd[vtmp.rc];
	for(; pi<pend; pi++)
		for(pj=pi+1; pj<pend; pj++)
			if(*pi<*pj) swapt(*pi,*pj,double);
	return vtmp;
}

//***************************  class CIIR  *********************************/
CIIR::CIIR(void)
{
}

CIIR::CIIR(double *b0, double *a0, int n0)
{
	psinsassert(n0<IIRnMax);
	for(int i=0; i<n0; i++)  { b[i]=b0[i]/a0[0]; a[i]=a0[i]; x[i]=y[i]=0.0; }
	n = n0;
}

double CIIR::Update(double x0)
{
//	a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
//                        - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
	double y0 = 0.0;
	for(int i=n-1; i>0; i--)
	{
		x[i] = x[i-1]; y[i] = y[i-1];
		y0 += b[i]*x[i] - a[i]*y[i];
	}
	x[0] = x0;
	y0 += b[0]*x0;
	y[0] = y0;
	return y0;
}

//***************************  class CV3IIR  *********************************/
CIIRV3::CIIRV3(void)
{
}

CIIRV3::CIIRV3(double *b0, double *a0, int n0, double *b1, double *a1, int n1, 
			   double *b2, double *a2, int n2)
{
	iir0 = CIIR(b0, a0, n0);
	if(n1==0)	iir1 = iir0;	// iir1 the same as iir0
	else		iir1 = CIIR(b1, a1, n1);
	if(n2==0)	iir2 = iir0;	// iir2 the same as iir0
	else		iir2 = CIIR(b2, a2, n2);
}

CVect3 CIIRV3::Update(const CVect3 &x)
{
	return y = CVect3(iir0.Update(x.i), iir1.Update(x.j), iir2.Update(x.k));
}

//***************************  class CRAvar  *********************************/
CRAvar::CRAvar()
{
}

CRAvar::CRAvar(int nR00, int maxCount0)
{
	psinsassert(nR00<RAMAX);
	this->nR0 = nR00;
	for(int i=0; i<RAMAX; i++)  { Rmaxcount[i]=maxCount0, tau[i]=INF; }
}

void CRAvar::set(double r00, double tau0, double rmax, double rmin, int i)
{
	this->R0[i] = r00*r00;
	this->tau[i] = tau0;
	this->r0[i] = 0.0;  Rmaxflag[i] = Rmaxcount[i];
	this->Rmax[i] = abs(rmax)<EPS ? 100.0*this->R0[i] : rmax*rmax;
	this->Rmin[i] = abs(rmin)<EPS ?  0.01*this->R0[i] : rmin*rmin;
}

void CRAvar::set(const CVect3 &r00, const CVect3 &tau0, const CVect3 &rmax, const CVect3 &rmin)
{
	const double *pr0=&r00.i, *ptau=&tau0.i, *prmax=&rmax.i, *prmin=&rmin.i;
	for(int i=0; i<3; i++,pr0++,ptau++,prmax++,prmin++)
		set(*pr0, *ptau, *prmax, *prmin, i);
}

void CRAvar::set(const CVect &r00, const CVect &tau0, const CVect &rmax, const CVect &rmin)
{
	const double *pr0=r00.dd, *ptau=tau0.dd, *prmax=rmax.dd, *prmin=rmin.dd;
	for(int i=0; i<nR0; i++,pr0++,ptau++,prmax++,prmin++)
		set(*pr0, *ptau, *prmax, *prmin, i);
}

void CRAvar::Update(double r, double ts0, int i)
{
	if(tau[i]>INFp5) return;
	double tstau = ts0>tau[i] ? 1.0 : ts0/tau[i];
	double dr2=r-r0[i]; dr2=dr2*dr2; r0[i]=r;
	if(dr2>R0[i]) R0[i]=dr2; else R0[i]=(1.0-tstau)*R0[i]+tstau*dr2;
	if(R0[i]<Rmin[i]) R0[i]=Rmin[i];
	if(R0[i]>Rmax[i]) {R0[i]=Rmax[i];Rmaxflag[i]=Rmaxcount[i];} else {Rmaxflag[i]-=Rmaxflag[i]>0;}
}

void CRAvar::Update(const CVect3 &r, double ts0)
{
	const double *pr=&r.i;
	for(int i=0; i<3; i++,pr++)
		Update(*pr, ts0, i);
}

void CRAvar::Update(const CVect &r, double ts0)
{
	const double *pr=r.dd;
	for(int i=0; i<nR0; i++,pr++)
		Update(*pr, ts0, i);
}

double CRAvar::operator()(int k) const
{
	return Rmaxflag[k] ? INF : sqrt(R0[k]);
}

//***************************  class CVAR  ***********************************/
CVAR::CVAR(int imax0, double data0)
{
	imax = mmin(imax0, VARMAX);
	for(ipush=0; ipush<imax; ipush++)	array[ipush] = data0;
	ipush = 0;
	mean = data0; var = 0.0;
}

double CVAR::Update(double data, BOOL isvar)
{
	array[ipush] = data;
	if(++ipush==imax) ipush=0;
	double *p0, *p1;
	for(mean=0.0,p0=&array[0],p1=&array[imax]; p0<p1; p0++)  mean+=*p0;
	mean /= imax;
	if (isvar)
	{
		for (var = 0.0, p0 = &array[0], p1 = &array[imax]; p0 < p1; p0++)  { double vi = *p0 - mean; var += vi*vi; }
		var /= imax - 1;
	}
	return var;
}

//***************************  class CVARn  ***********************************/
CVARn::CVARn(void)
{
	pData = NULL;
}

CVARn::CVARn(int row0, int clm0)
{
	row = row0, clm = clm0;
	pData = new double*[clm];
	pData[0] = new double[row*clm+5*clm];
	for(int i = 1; i < clm; i++) {
		pData[i] = pData[i - 1] + row;
	}
	pd = pData[clm-1]+row;  Sx = pd + clm;  Sx2 = Sx + clm;  mx = Sx2 + clm;  stdx = mx + clm;
	stdsf = sqrt((row - 1) / (row - 2.0));
	Reset();
}

CVARn::~CVARn(void)
{
	Deletep(pData);
}

void CVARn::Reset(void)
{
	idxpush = rowcnt = 0;
	memset(pData[0], 0, row*clm*sizeof(double));
	memset(Sx, 0, 5*clm*sizeof(double));
}

BOOL CVARn::Update(const double *pf)
{
	if (!pData[0])  return FALSE;
	if (++rowcnt > row) rowcnt = row;
	int idxnext = (idxpush >= row - 1) ? 0 : idxpush + 1;
	for(int i = 0; i < clm; i++)
	{
		double f=*pf++;  if(f>1e5) f=1e5; else if(f<-1e5) f=-1e5;
		pData[i][idxpush] = f;
		Sx[i] += f - pData[i][idxnext];
		Sx2[i] += f*f - pData[i][idxnext] * pData[i][idxnext];
		mx[i] = Sx[i] / rowcnt;
		stdx[i] = sqrt(Sx2[i] / rowcnt - mx[i] * mx[i]) * stdsf;   // Dx = E(x^2) - (Ex)^2
	}
	if (++idxpush == row) {
		idxpush = 0;
	}
	return idxpush == 0;
}

BOOL CVARn::Update(double f, ...)
{
	va_list vl;
	va_start(vl, f);
	for(int i = 0; i < clm; i++)
	{
		pd[i] = f;
		f = va_arg(vl, double);
	}
	va_end(vl);
	return Update(pd);
}

//***************************  class CContLarge  *********************************/
CContLarge::CContLarge(void)
{
}

CContLarge::CContLarge(double large0, double dt0, int cnt0)
{
	large=large0; dt=dt0; cnt=cnt0; cnti=0; t0=-INF;
}
    
BOOL CContLarge::Update(double val, double t1)
{
	if(t1<INFp5) {
		if(t1-t0>dt) cnti=0;  // restart
		t0=t1;
	}
	if(val>large||val<-large) {
		if(++cnti>cnt)	{ cnti=0; return TRUE; }
	}
	return FALSE;
};

//***************************  class CAbnomalCnt  *********************************/
CAbnomalCnt::CAbnomalCnt(void)
{
}

CAbnomalCnt::CAbnomalCnt(int cntMax0, double tlast0, double valMax0, double valMin0)
{
	cntMax = cntMax0;  tlast = tlast0;  valMax = valMax0;  valMin = valMin0<-INFp5 ? -valMax : valMin0;
	cnt = abnFlag = 0;  t0 = 0.0;
}

BOOL CAbnomalCnt::Update(double val, double t)
{
	if(abnFlag && (++cnt>cntMax||(t-t0)>tlast)) abnFlag = FALSE;
	if(val<valMin||val>valMax)	{
		abnFlag = TRUE;  cnt = 0;  t0 = t;
	}
	return abnFlag;
}

//***************************  class CWzhold  *********************************/
CWzhold::CWzhold(void)
{
}

void CWzhold::Init(double maxw0, double ts0, double T0, int cntNP0)
{
	ts = ts0;  T = T0;	maxw = maxw0;	cntNP = cntNP0;
	meanw = meanwpre = meanwprepre = tstop = 0.0;
	Reset();
}
	
void CWzhold::Reset(void)
{
	t = val = 0.0;
	cntNi = cntPi = big = 0;
}
	
int CWzhold::Update(double wz)   // wz in rad/s
{
	retres=0;
	if((tstop+=ts)<3.0) {
		retres=-6;
	}
	else {
		if(wz<-maxw) {
			if(++cntNi>cntNP) { meanw=meanwpre=meanwprepre=tstop=0.0; Reset(); retres=-3; }
			else retres=-1;
			big++;
		}
		else if(wz>maxw) {
			if(++cntPi>cntNP) { meanw=meanwpre=meanwprepre=tstop=0.0; Reset(); retres=-4; }
			else retres=-2;
			big++;
		}
		else {
			cntNi = cntPi = 0;
			big--;
		}
		if(big<0) {
			big=0;
		}
		else if(big>10) {
			meanw=meanwpre=meanwprepre=0.0; Reset(); retres=-5;
		}
		if(retres>=-2) {
			val += wz;
			t += ts;
			if(t>T) {
				meanwprepre = meanwpre;
				meanwpre = meanw;
				meanw = val*ts/T;  // meanw in rad/s
				Reset();
				if(meanwpre<-EPS||meanwpre>EPS) {
					if(meanwprepre<-EPS||meanwprepre>EPS) retres=3;
					else retres = 2;
				}
				else {
					retres = 1;
				}
			}
		}
	}
	return retres;
};

//***************************  class CMaxMin  *********************************/
CMaxMin::CMaxMin(int cnt00, int pre00, float f0)
{
	Init(cnt00, pre00, f0);
}

void CMaxMin::Init(int cnt00, int pre00, float f0)
{
	max0=f0, min0=-f0, maxpre0=f0, minpre0=-f0;
	maxCur=f0, minCur=-f0, maxpreCur=f0, minpreCur=-f0;
	maxRes=f0, minRes=-f0, maxpreRes=f0, minpreRes=-f0, diffRes=diffpreRes=2*f0, meanRes=sumRes=maxabsRes=0.0;
	cntCur=cnt0=cnt00;
	cntpreCur = (pre00<=0||pre00>=cnt00) ? cnt00/2 : cnt0-pre00;
	flag = 0;
}

void CMaxMin::Restart(void)
{
	Init(cnt0);
}

int CMaxMin::Update(float f)
{
	flag=0;
	if(maxCur<f) maxCur=f; else if(minCur>f) minCur=f;
	if(maxpreCur<f) maxpreCur=f; else if(minpreCur>f) minpreCur=f;
	sumRes += f;
	if(--cntCur<=0) {
		maxRes=maxCur; minRes=minCur; maxCur=minCur=f; diffRes=maxRes-minRes; meanRes=sumRes/cnt0; cntCur=cnt0; flag=1;
		maxabsRes = minRes<-maxRes ? -minRes : maxRes;
		sumRes = 0.0;
	}
	if(--cntpreCur<=0) {
		maxpreRes=maxpreCur; minpreRes=minpreCur; maxpreCur=minpreCur=f; diffpreRes=maxpreRes-minpreRes; cntpreCur=cnt0; flag=-1;
	}
	return flag;
}

//***************************  class CMaxMinn  *********************************/
CMaxMinn::CMaxMinn(int n0, int cnt00, int pre00, float f0)
{
	Init(n0, cnt00, pre00, f0);
}

void CMaxMinn::Init(int n0, int cnt00, int pre00, float f0)
{
	n = n0;
	CMaxMin mm0(cnt00, pre00, f0);
	for(int i=0; i<n; i++)	mm[i] = mm0;
	flag = 0;
}

void CMaxMinn::Restart(void)
{
	for(int i=0; i<n; i++) mm[i].Restart();
	flag = 0;
}

int CMaxMinn::Update(float f, ...)
{
	CMaxMin *pm=&mm[0];
	va_list vl;
	va_start(vl, f);
	for(int i=0; i<n; i++,pm++)
	{
		pm->Update(f);
		f = (float)va_arg(vl, double);
	}
	va_end(vl);
	return flag = mm[0].flag;
}

int CMaxMinn::Update(const CVect3 &v1)
{
	mm[0].Update((float)v1.i); mm[1].Update((float)v1.j); mm[2].Update((float)v1.k);
	return flag = mm[0].flag;
}

int CMaxMinn::Update(const CVect3 &v1, const CVect3 &v2)
{
	mm[0].Update((float)v1.i); mm[1].Update((float)v1.j); mm[2].Update((float)v1.k);
	mm[3].Update((float)v2.i); mm[4].Update((float)v2.j); mm[5].Update((float)v2.k);
	return flag = mm[0].flag;
}

int CMaxMinn::Update(const CVect3 &v1, const CVect3 &v2, const CVect3 &v3)
{
	mm[0].Update((float)v1.i); mm[1].Update((float)v1.j); mm[2].Update((float)v1.k);
	mm[3].Update((float)v2.i); mm[4].Update((float)v2.j); mm[5].Update((float)v2.k);
	mm[6].Update((float)v3.i); mm[7].Update((float)v3.j); mm[8].Update((float)v3.k);
	return flag = mm[0].flag;
}

int CMaxMinn::Update(const CVect3 &v1, const CVect3 &v2, const CVect3 &v3, const CVect3 &v4)
{
	mm[0].Update((float)v1.i); mm[1].Update((float)v1.j); mm[2].Update((float)v1.k);
	mm[3].Update((float)v2.i); mm[4].Update((float)v2.j); mm[5].Update((float)v2.k);
	mm[6].Update((float)v3.i); mm[7].Update((float)v3.j); mm[8].Update((float)v3.k);
	mm[9].Update((float)v4.i); mm[10].Update((float)v4.j); mm[11].Update((float)v4.k);
	return flag = mm[0].flag;
}

int CMaxMinn::Update(const CVect3 &v1, const CVect3 &v2, const CVect3 &v3, const CVect3 &v4, const CVect3 &v5)
{
	mm[0].Update((float)v1.i); mm[1].Update((float)v1.j); mm[2].Update((float)v1.k);
	mm[3].Update((float)v2.i); mm[4].Update((float)v2.j); mm[5].Update((float)v2.k);
	mm[6].Update((float)v3.i); mm[7].Update((float)v3.j); mm[8].Update((float)v3.k);
	mm[9].Update((float)v4.i); mm[10].Update((float)v4.j); mm[11].Update((float)v4.k);
	mm[12].Update((float)v5.i); mm[13].Update((float)v5.j); mm[14].Update((float)v5.k);
	return flag = mm[0].flag;
}

float CMaxMinn::ResFloat(int i, int minmeanmaxFlag)
{
	     if(minmeanmaxFlag==-1) return mm[i].minRes;    // min
	else if(minmeanmaxFlag==0 ) return mm[i].meanRes;   // mean
	else if(minmeanmaxFlag==1 ) return mm[i].maxRes;    // max
	else if(minmeanmaxFlag==2 ) return mm[i].maxabsRes; // maxabs
	else return 0.0f;
}

CVect3 CMaxMinn::ResVect3(int i, int minmeanmaxFlag)
{
	CMaxMin *pi0=&mm[i], *pi1=pi0+1, *pi2=pi0+2;
	     if(minmeanmaxFlag==-1) return CVect3(pi0->minRes,    pi1->minRes,    pi2->minRes);    // min
	else if(minmeanmaxFlag==0 ) return CVect3(pi0->meanRes,   pi1->meanRes,   pi2->meanRes);   // mean
	else if(minmeanmaxFlag==1 ) return CVect3(pi0->maxRes,    pi1->maxRes,    pi2->maxRes);    // max
	else if(minmeanmaxFlag==2 ) return CVect3(pi0->maxabsRes, pi1->maxabsRes, pi2->maxabsRes); // maxabs
	else if(minmeanmaxFlag==3 ) return CVect3(pi0->diffRes,   pi1->diffRes,   pi2->diffRes);   // diffRes
	else return O31;
}

//***************************  class CKalman  *********************************/
CKalman::CKalman(void)
{
}

CKalman::CKalman(int nq0, int nr0)
{
//	psinsassert(nq0<=MMD&&nr0<=MMD);
	if(!(nq0<=MMD&&nr0<=MMD)) { /*printf("\tMMD too small!\n");*/ exit(0); }
	Init(nq0, nr0);
}

void CKalman::Init(int nq0, int nr0)
{
	kftk = 0.0;
	nq = nq0; nr = nr0;
	Ft = Pk = CMat(nq,nq,0.0);
	Hk = CMat(nr,nq,0.0);  Fading = CMat(nr,nq,1.0); zfdafa = 0.1f;
	Qt = Pmin = Xk = CVect(nq,0.0);  Xmax = Pmax = CVect(nq,INF);  Pset = CVect(nq,-INF);
	Zk = Zk_1 = CVect(nr,0.0);  Rt = Rt0 = CVect(nr,INF); Rset = CVect(nr,-INF); rts = CVect(nr,1.0);  Zfd = CVect(nr,0.0); Zfd0 = Zmm0 = Zmax = CVect(nr,INF);
	innoPre = CVect(nr,0.0); innoDiffMax = CVect(nr,INF);
	RtTau = Rmax = CVect(nr,INF); measstop = measlost = Rmin = Rb = Rstop = CVect(nr,0.0); Rbeta = CVect(nr,1.0);
	SetRmaxcount(5);
	innoMax = CVect(nr,INF);
	SetInnoMaxcount(5);
	FBTau = FBMax = FBOne = FBOne1 = CVect(nq,INF); FBXk = FBTotal = CVect(nq,0.0);
	kfcount = measflag = measflaglog = 0;  SetMeasMask(nr0,3);
	Zmm.Init(nr0, 10);
}

void CKalman::SetRmmbt(double rmin, double rmax, double b, double tau)
{
	if(tau<INF/2)  RtTau = tau;
	if(b<INF/2)  Rb = b;
	if(rmax<INF/2)  Rmax = Rt * (rmax*rmax);
	Rmin = Rt * (rmin*rmin);
}

void CKalman::SetRmaxcount(int cnt)
{
	for(int i=0; i<nr; i++) { Rmaxcount[i]=0, Rmaxcount0[i]=cnt; }
}

void CKalman::SetInnoMaxcount(int cnt)
{
	for(int i=0; i<nr; i++) { innoMaxcount[i]=0, innoMaxcount0=cnt; }
}

void CKalman::SetZmm(int zi, int pi, double zmm0, int cnt)
{
	if(cnt>0) Zmm.mm[zi].Init(cnt);
	Zmmpk[zi] = pi;  Zmm0.dd[zi] = zmm0; 
}

void CKalman::SetZmmVn(const CVect3 &zmm0, int cnt)
{
	SetZmm(0, 3, zmm0.i, cnt);	SetZmm(1, 4, zmm0.j, cnt);	SetZmm(2, 5, zmm0.k, cnt);
}

void CKalman::SetZmmPos(const CVect3 &zmm0, int cnt)
{
	SetZmm(3, 6, zmm0.i, cnt);	SetZmm(4, 7, zmm0.j, cnt);	SetZmm(5, 8, zmm0.k, cnt);
}

void CKalman::TimeUpdate(double kfts0, int fback)
{
	CMat Fk;
	kftk += kfts0;  kfcount++;
	SetFt(nq);
	Fk = ++(Ft*kfts0);  // Fk = I+Ft*ts
	Xk = Fk * Xk;
	Pk = Fk*Pk*(~Fk);  Pk += Qt*kfts0;
	if(fback)  Feedback(nq, kfts0);
	for(int i=0; i<nr; i++) {
		measlost.dd[i] += kfts0;
		if(abs(measstop.dd[i])<EPS) measstop.dd[i] -= kfts0;
	}
}

void CKalman::SetMeasFlag(unsigned int flag, int type)
{
	if(type==1)
		measflag = (flag==0) ? 0 : (measflag|flag);  // add the flag bits to 1
	else if(type==0)
		measflag &= ~flag;  // set the flag bits to 0
}

void CKalman::SetMeasMask(unsigned int mask, int type)
{
	int m;
	if(type==1) measmask = mask;		// set mask 1
	else if(type==0) measmask &= ~mask;	// delete mask to 0
	else if(type==2) measmask |= mask;	// add mask 1
	else if(type==3) {					// set mask-LSB 1
		for(m=0; mask>0; mask--)  m |= 1<<(mask-1);
		SetMeasMask(m);
	}
}

void CKalman::SetMeasStop(unsigned int meas, double stop)
{
	for(int i=0; i<measstop.rc; i++)  {	// assert(rc<32)
		if( meas&(0x01<<i) && (measstop.dd[i]<stop||stop<=0.0) ) measstop.dd[i]=stop;
	}
//	measstop.SetBit(meas, stop);
}

void CKalman::SetRadptStop(unsigned int meas, double stop)
{
	Rstop.SetBit(meas, stop);
}

int CKalman::MeasUpdate(double fading)
{
	CVect Pxz, Kk, Hi;
	SetMeas();
	for(int i=0; i<nr; i++)
	{
		if(((measflag&measmask)&(0x01<<i)) && measstop.dd[i]>EPS)
		{
			Hi = Hk.GetRow(i);
			Pxz = Pk*(~Hi);
			double Pz0 = (Hi*Pxz)(0,0), r=Zk(i)-(Hi*Xk)(0,0);
			if(Rb.dd[i]>EPS)
				RAdaptive(i, r, Pz0);
			if(Zfd.dd[i]<INFp5)
				RPkFading(i);
			double Pzz = Pz0+Rt.dd[i]/rts.dd[i];
			Kk = Pxz*(1.0/Pzz);
			Xk += Kk*r;
			Pk -= Kk*(~Pxz);
			measlost.dd[i] = 0.0;
		}
	}
	if(fading>1.0) Pk *= fading;
	XPConstrain();
	symmetry(Pk);
	int measres = measflag&measmask;
	measflaglog |= measres;
	SetMeasFlag(0);
	return measres;
}

void MeasUD(double U[], double D[], const double H[], double R, double K[], int n)  // Ref: my book P291
{
	int i, j, s;
	double f[MMD], g[MMD], DUHR[MMD], afa=R;
	for(i=0; i<n; i++) {
		for(f[i]=H[i],j=0; j<i; j++) f[i] += H[j]*U[j*n+i];  // U: upper triangle 
		g[i] = D[i]*f[i];
		afa += f[i]*g[i];
	}
	for(j=n-1; j>=0; j--) {
		double afa0=afa-f[j]*g[j], lambda=-f[j]/afa0;
		DUHR[j] = H[j];
		D[j] *= afa0/afa;  afa=afa0;
		for(i=j-1; i>=0; i--) {
			U[i*n+j] += lambda*g[i]; 
			for(s=i+1; s<j; s++) U[i*n+j] += lambda*U[i*n+s]*g[s];
			DUHR[j] += H[i]*U[i*n+j];  // H*U = U^T*H^T
		}
		DUHR[j] *= D[j]/R;
	}
	for(i=0; i<n; i++) {
		for(K[i]=DUHR[i],j=i+1; j<n; j++) K[i] += U[i*n+j]*DUHR[j];  // K=U*D*U^T*H^T/R
	}
}

int CKalman::RAdaptive(int i, double r, double Pr)
{
	double rr=r*r-Pr;
	if(Rb.dd[i]>1.0)	{  // s^2=Rb.dd[i], for heavy-tailed noise
		Rt.dd[i] = rr/Rb.dd[i];  // rho^2/s^2;
		if(Rt.dd[i]<Rmin.dd[i]) Rt.dd[i]=Rmin.dd[i];
		return 1; //adptOK=1;
	}
	if(rr<Rmin.dd[i])	rr = Rmin.dd[i];
	if(rr>Rmax.dd[i])	{ Rt.dd[i]=Rmax.dd[i]; Rmaxcount[i]++; }  
	else				{ Rt.dd[i]=(1.0-Rbeta.dd[i])*Rt.dd[i]+Rbeta.dd[i]*rr; Rmaxcount[i]=0; }
	Rbeta.dd[i] = Rbeta.dd[i]/(Rbeta.dd[i]+Rb.dd[i]);   // beta = beta / (beta+b)
	int adptOK = (Rmaxcount[i]==0||Rmaxcount[i]>Rmaxcount0[i]) ? 1: 0;
	return adptOK;
}

void CKalman::RPkFading(int i)
{
	Zfd.dd[i] = Zfd.dd[i]*((double)1.0f-zfdafa) + Zk.dd[i]*zfdafa;
	if(Zfd.dd[i]>Zfd0.dd[i] || Zfd.dd[i]<-Zfd0.dd[i])
		DVMDVafa(Fading.GetRow(i), Pk);
}

void CKalman::ZmmPkSet(int i)
{
	CMaxMin *pmm=&Zmm.mm[i];
	pmm->Update((float)Zk.dd[i]);
	if(pmm->flag) {  // if conitnously big abs(Zk[i]), enlarge Pk[i,i]
		if( (double)pmm->minRes>Zmm0.dd[i] || (double)pmm->maxRes<-Zmm0.dd[i] )
			Pset.dd[Zmmpk[i]] = pmm->meanRes*pmm->meanRes;
	}
}

void CKalman::XPConstrain(void)
{
	int i=0, nq1=nq+1;
	for(double *px=Xk.dd,*pxmax=Xmax.dd,*p=Pk.dd,*pmin=Pmin.dd,*pminEnd=&Pmin.dd[nq],*pmax=Pmax.dd,*pset=Pset.dd;
		pmin<pminEnd; px++,pxmax++,p+=nq1,pmin++,pmax++,pset++)
	{
		if(*px>*pxmax)		// Xk constrain
		{
			*px = *pxmax;
		}
		else if(*px<-*pxmax)
		{
			*px = -*pxmax;
		}
		if(*p<*pmin)		// Pk constrain
		{
			*p = *pmin;
		}
		else if(*p>*pmax)
		{
			double sqf=sqrt(*pmax/(*p))*0.9;
			for(double *prow=&Pk.dd[i*Pk.clm],*prowEnd=prow+nq,*pclm=&Pk.dd[i]; prow<prowEnd; prow++,pclm+=nq)
			{
				*prow *= sqf;
				*pclm = *prow;
			}
			Pk.dd[i*Pk.clm+i] *= sqf;  //20200303
			break;
		}
		if(*pset>0.0)	// Pk set
		{
			if(*p<*pset)
			{
				*p = *pset;
			}
			else if(*p>*pset)
			{
				double sqf=sqrt(*pset/(*p));
				for(double *prow=&Pk.dd[i*Pk.clm],*prowEnd=prow+nq,*pclm=&Pk.dd[i]; prow<prowEnd; prow++,pclm+=nq)
				{
					*prow *= sqf;
					*pclm *= sqf;
				}
			}
			*pset = -1.0;
		}
		i++;
	}
}

void CKalman::PmaxPminCheck(void)
{
	for(double *p=Pk.dd,*pmax=Pmax.dd,*pmin=Pmin.dd,*pminEnd=&Pmin.dd[nq]; pmin<pminEnd; p+=nq+1,pmax++,pmin++)
	{
		if(*p>*pmax) *pmax = *p*10.0;
		if(*p<EPS)	 *pmin = 0.0;  else if(*p<*pmin)  *pmin = *p/2.0;
	}
}

void CKalman::Feedback(int nnq, double fbts)
{
	double *pTau=FBTau.dd, *pTotal=FBTotal.dd, *pMax=FBMax.dd, *pOne=FBOne.dd, *pOne1=FBOne1.dd, *pXk=FBXk.dd, *p=Xk.dd;
	for(int i=0; i<nq; i++, pTau++,pTotal++,pMax++,pOne++,pOne1++,pXk++,p++)
	{
		if(*pTau<INFp5)
		{
			if(*p>*pOne1 || *p<-*pOne1) {  // max feedback/one step
				*pXk=*p;
			}
			else {
				double afa = fbts<*pTau ? fbts/(*pTau) : 1.0;
				*pXk = *p*afa;
				if(*pXk>*pOne) *pXk=*pOne; else if(*pXk<-*pOne) *pXk=-*pOne;  // min feedback/one step
			}
			if(*pMax<INFp5)
			{
				if(*pTotal+*pXk>*pMax)			*pXk = *pMax-*pTotal;
				else if(*pTotal+*pXk<-*pMax)	*pXk = -*pMax-*pTotal;
			}
			*p -= *pXk;
			*pTotal += *pXk;
		}
		else
		{
			*pXk = 0.0;
		}
	}
}

void CKalman::FeedbackAll(void)
{
	CVect oldFB=FBTau, oldFBOne=FBOne, oldFBMax=FBMax;
	FBTau = 0.0;  FBOne = INF;  FBMax = INF;
	Feedback(nq, 1.0);
	FBTau = oldFB;  FBOne = oldFBOne;  FBMax = oldFBMax;
	measflag = 0;
}

void CKalman::RtFading(int i, double fdts)
{
	double Taui=RtTau.dd[i], Rti=Rt.dd[i], Rmaxi=Rmax.dd[i];
	if(measlost.dd[i]>3.0 && Taui<INFp5 && Rb.dd[i]>0.0 && Rti<Rmaxi)
	{
		double afa = fdts<Taui ? fdts/Taui : 1.0;
		Rti += 2*sqrt(Rmaxi*Rti)*afa;
		Rt.dd[i] = Rti;
	}
}

void fusion(double *x1, double *p1, const double *x2, const double *p2, int n, double *xf, double *pf)
{
	if(xf==NULL) { xf=x1, pf=p1; }
	double *x10=nullptr, *xf0=nullptr;
	CVect3 att1;
	if(n<100) {  // n<100 for att(1:3), n>100 for no_att(1:3)
		x10 = x1, xf0 = xf;
		att1 = *(CVect3*)x1; 
		*(CVect3*)x2 = qq2phi(a2qua(*(CVect3*)x2),a2qua(att1));
		*(CVect3*)x1 = O31;
	}
	int j;
	for(j=(n>100)?100:0; j<n; j++,x1++,p1++,x2++,p2++,xf++,pf++)
	{
		double p1p2 = *p1+*p2;
		*xf = (*p1**x2 + *p2**x1)/p1p2; 
		*pf = *p1**p2/p1p2;
	}
	if(j<100) {
		*(CVect3*)xf0 = q2att(a2qua(att1)+*(CVect3*)xf0);
		if(xf0!=x10) *(CVect3*)x10 = att1; 
	}
}

void fusion(CVect3 &x1, CVect3 &p1, const CVect3 x2, const CVect3 p2)
{
	fusion(&x1.i, &p1.i, &x2.i, &p2.i, 103);
}

void fusion(CVect3 &x1, CVect3 &p1, const CVect3 x2, const CVect3 p2,
			CVect3 &xf, CVect3 &pf)
{
	fusion(&x1.i, &p1.i, &x2.i, &p2.i, 103, &xf.i, &pf.i);
}


//***************************  class CSINSTDKF  *********************************/
CSINSTDKF::CSINSTDKF(void)
{
}

CSINSTDKF::CSINSTDKF(int nq0, int nr0)
{
	CKalman::Init(nq0, nr0);
	Kkp2y = Kkv2y = 1.0;
}

void CSINSTDKF::Init(const CSINS &sins0)
{
	sins = sins0;  kftk = sins.tk;
	Fk = eye(nq);  Pk1 = CMat(nq,nq, 0.0);
	Pxz = Qk = Kk = tmeas = CVect(nq, 0.0);
	meantdts = 1.0; tdts = 0.0;
	maxStep = 2*(nq+nr)+3;
	TDReset();
	curOutStep = 0; maxOutStep = 1;
	timcnt0 = timcnt1 = 0, timcntmax = 100;  burden = 0.0;
	cststt = nq;  for(int k=0; k<nq; k++) hi1[k]=-1;
}

void CSINSTDKF::TDReset(void)
{
	iter = -2;
	ifn = 0;	meanfn = O31;
	SetMeasFlag(0);
}

double CSINSTDKF::Innovationi(int row)
{
	double hx=0.0;
	for(double *ph=&Hk.dd[row*nq], *px=&Xk.dd[0], *pxEnd=&Xk.dd[Xk.rc]; px<pxEnd; ph++,px++)  { hx += (*ph)*(*px); }
	return (Zk.dd[row]-hx);
}

int CSINSTDKF::TDUpdate(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, int nStep)
{
	sins.Update(pwm, pvm, nSamples, ts);
	if(++curOutStep>=maxOutStep) { RTOutput(), curOutStep=0; }
	Feedback(nq, sins.nts);
	for(int j=0; j<nr; j++) {
		measlost.dd[j] += sins.nts;
		if(Rstop.dd[j]>0.0) Rstop.dd[j] -= sins.nts;
		if(measstop.dd[j]>0.0) measstop.dd[j] -= sins.nts;
	}

	measRes = 0;

	if(nStep<=0||nStep>=maxStep) { nStep=maxStep; }
	tdStep = nStep;

	tdts += sins.nts; kftk = sins.tk;  kfcount++;
//	meanfn = meanfn+sins.fn; ifn++;
	VADDE(meanfn, sins.fn); ifn++;
	for(int i=0; i<nStep; i++)
	{
		if(iter==-2)			// -2: set measurements
		{
			if(ifn==0)	break;
			CVect3 fn=sins.fn, an=sins.an;
			sins.fn = meanfn*(1.0/ifn); meanfn = O31; ifn = 0;
			sins.an = sins.fn+sins.eth.gcc; 
			SetFt(nq);
			SetMeas(); SetHk(nq); sins.fn = fn; sins.an = an;
		}
		else if(iter==-1)			// -1: discrete
		{
//			Fk = ++(Ft*tdts); // Fk = I+Ft*ts
			double *pFk,*pFt,*pEnd;
			for(pFk=Fk.dd,pFt=Ft.dd,pEnd=&Fk.dd[cststt*Fk.clm]; pFk<pEnd; pFk++,pFt++)  *pFk=*pFt*tdts;
			for(pFk=Fk.dd; pFk<pEnd; pFk+=Fk.clm+1)  *pFk+=1.0;
//			Xk = Fk*Xk;
			pFk=Fk.dd, pEnd=&Xk.dd[Xk.rc];  int jj;
			for(jj=0; jj<cststt; jj++)  {
				double f=0.0;
				for(double *pX=Xk.dd; pX<pEnd; pFk++,pX++)  f += *pFk * *pX;
				Pxz.dd[jj] = f;  // Pxz just for store data
			}
			for(jj=0; jj<cststt; jj++)  Xk.dd[jj] = Pxz.dd[jj];
#ifndef PSINS_FAST_CALCULATION
			Qk = Qt*tdts;
#else
			mul(Qk, Qt, tdts);
#endif
//			RtFading(tdts);
			meantdts = tdts; tdts = 0.0;
		}
		else if(iter<nq)		// 0 -> (nq-1): Fk*Pk
		{
			int row=iter;
			RowMul(Pk1, Fk, Pk, row, cststt);
		}
		else if(iter<2*nq)		// nq -> (2*nq-1): Fk*Pk*Fk+Qk
		{
			int row=iter-nq;
			RowMulT(Pk, Pk1, Fk, row, cststt);
			Pk.dd[nq*row+row] += Qk.dd[row];
//			if(row==nq-1) {	Pk += Qk; }
		}
		else if(iter<2*(nq+nr))	// (2*nq) -> (2*(nq+nr)-1): sequential measurement updating
		{
			int row=(iter-2*Ft.row)/2;
			int flag = (measflag&measmask)&(0x01<<row);
			if(flag)
			{
//				if((iter-2*Ft.row)%2==0)
				if(iter%2==0)
				{
					Hk.GetRow(Hi, row);  int hi=hi1[row];
					if(hi>=0) {
						Pk.GetClm(Pxz, hi);
						Pz0 = Pxz.dd[hi];
						innovation = Zk.dd[row]-Xk.dd[hi];  Zk_1.dd[row] = Zk.dd[row];
					}
					else {
#ifndef PSINS_FAST_CALCULATION
						Pxz = Pk*(~Hi);
#else
						mul(Pxz, Pk, Hi);
#endif
						Pz0 = dot(Hi,Pxz);
						innovation = Zk.dd[row]-dot(Hi,Xk);  Zk_1.dd[row] = Zk.dd[row];
					}
					double innoDiff = innovation - innoPre.dd[row];  innoPre.dd[row] = innovation;
					int *pinnoMaxcount=&innoMaxcount[row];
					if(*pinnoMaxcount>0) (*pinnoMaxcount)--;
					if(innovation<-innoMax.dd[row] || innovation>innoMax.dd[row])
					{
						if(*pinnoMaxcount<2*innoMaxcount0) *pinnoMaxcount += 2;
					}
					if( (*pinnoMaxcount==0||*pinnoMaxcount>innoMaxcount0) 
						&& (innovation>-Zmax.dd[row]&&innovation<Zmax.dd[row])
						&& (innoDiff>-innoDiffMax.dd[row]&&innoDiff<innoDiffMax.dd[row]) )
					{
						adptOKi = 1;
						if(Rb.dd[row]>0.0 && Rstop.dd[row]<=0.0) {
							adptOKi=RAdaptive(row, innovation, Pz0);
						}
						// if(row==5)
						// 	int a=1;
						if(Rset.dd[row]>0.0) {
							Rt.dd[row]=Rset.dd[row];  Rset.dd[row]=-1.0;  adptOKi=1;
						}
						double Pzz = Pz0 + Rt.dd[row]/rts.dd[row];
#ifndef PSINS_FAST_CALCULATION
						Kk = Pxz*(1.0/Pzz);
#else
						mul(Kk, Pxz, (1.0/Pzz));
#endif
					}
					else
					{
						adptOKi = 0;
					}
				}
				else
				{
					measflag ^= flag;
					if(adptOKi && measstop.dd[row]<EPS)
					{
						measRes |= flag;
#ifndef PSINS_FAST_CALCULATION
						Pk -= Kk*(~Pxz);
#else
						double *pPk=Pk.dd, *pKk=Kk.dd, *pKkEnd=&Kk.dd[Kk.row], *pPxzEnd=&Pxz.dd[Pxz.row];
						for(; pKk<pKkEnd; pKk++) {
							for(double *pPxz=Pxz.dd; pPxz<pPxzEnd; pPxz++)
								*pPk++ -= *pKk * *pPxz;
						}
#endif
						if(flag&030) Kk.dd[2]*=Kkp2y;  // disable pos2yaw
						else if(flag&003) Kk.dd[2]*=Kkv2y;
#ifndef PSINS_FAST_CALCULATION
						Xk += Kk*innovation;
#else
						for(double *pXk=Xk.dd,*pXkEnd=&Xk.dd[Xk.row],*pKk1=Kk.dd; pXk<pXkEnd; pXk++,pKk1++) *pXk += *pKk1 * innovation;
#endif
						measlost.dd[row] = 0.0;
					}
					if(Zfd0.dd[row]<INFp5)
					{
						RPkFading(row);
					}
					if(Zmm0.dd[row]<INFp5)
					{
						ZmmPkSet(row);
					}
				}
			}
			else
			{
				nStep++;
			}
			if(iter%2==0 && Rstop.dd[row]<EPS)
				RtFading(row, meantdts);
		}
		else if(iter==2*(nq+nr))	// 2*(nq+nr): Xk,Pk constrain & symmetry
		{
			XPConstrain();
			symmetry(Pk);
		}
		else if(iter>=2*(nq+nr)+1)	// 2*(nq+nr)+1: Miscellanous
		{
			Miscellanous();
			iter = -3;
		}
		iter++;
	}
	SecretAttitude();

	measflaglog |= measRes;
	return measRes;
}

void CSINSTDKF::MeasUpdate(const CVect &Hi0, double Ri, double Zi)
{
	if(iter>=0 && iter<2*nq) return;  // !
	Pxz = Pk*(~Hi0);
	Pz0 = (Hi0*Pxz)(0,0);
	innovation = Zi-(Hi0*Xk)(0,0);
	double Pzz = Pz0 + Ri;
	Kk = Pxz*(1.0/Pzz);
	Xk += Kk*innovation;
	Pk -= Kk*(~Pxz);
}

void CSINSTDKF::MarkovGyro(const CVect3 &tauG, const CVect3 &sRG, int stateeb)
{
	sins.SetTauGA(tauG, O31);
	*(CVect3*)&Qt.dd[stateeb] = MKQt(sRG, sins.tauGyro);
}

void CSINSTDKF::MarkovAcc(const CVect3 &tauA, const CVect3 &sRA, int statedb)
{
	sins.SetTauGA(O31, tauA);
	*(CVect3*)&Qt.dd[statedb] = MKQt(sRA, sins.tauAcc);
}

void CSINSTDKF::SetYaw(double yaw, int statephi, int statedvn)
{
	CQuat qnn = a2qua(0,0,diffYaw(yaw,sins.att.k));
	sins.qnb = qnn*sins.qnb;  sins.Cnb = q2mat(sins.qnb);  sins.Cbn = ~sins.Cnb;  sins.att = m2att(sins.Cnb);
	sins.vn = qnn*sins.vn;
	*(CVect3*)&Xk.dd[statephi] = qnn**(CVect3*)&Xk.dd[statephi];
	*(CVect3*)&Xk.dd[statedvn] = qnn**(CVect3*)&Xk.dd[statedvn];
	Pk = diag(diag(Pk));
/*	CMat3 Cnn=q2mat(qnn);
	CMat Cnn15(15,15,0.0);
	Cnn15.SetMat3(0,0, Cnn); Cnn15.SetMat3(3,3, Cnn);Cnn15.SetMat3(6,6, rcijk(Cnn,102)); 
	Cnn15.SetMat3(9,9, Cnn); Cnn15.SetMat3(12,12, Cnn);
	Pk = (Cnn15)*Pk*(~Cnn15);*/
	TDReset();
}

void CSINSTDKF::PSetVertCh(double sph, double spv, double spd)
{
	sph = sph*sph;  spv = spv*spv;  spd = spd*spd;
	if(sph>Pk.dd[nq*8+8])   Pset.dd[8]  = sph;
	if(spv>Pk.dd[nq*5+5])   Pset.dd[5]  = spv;
	if(spd>Pk.dd[nq*14+14]) Pset.dd[14] = spd;
}

double CSINSTDKF::SetCalcuBurden(unsigned int timcnt, int itype)
{
	double res=0.0;
	if(itype==-1) {
		timcntmax = timcnt;
	}
	else if(itype==0) {
		timcnt0 = timcnt;
	}
	else {
		timcnt1 = timcnt;
		if(timcntmax<timcnt1) timcntmax=timcnt1;
		if(timcnt1<timcnt0) timcnt1+=timcntmax;
		res = (double)(timcnt1-timcnt0)/timcntmax;
		if(res>0.99) res=0.99;
		burden = burden<res ? res : burden-0.01;
	}
	return (iter+100)+res;
}

//***************************  class CSINSGNSS  *********************************/
CSINSGNSS::CSINSGNSS(void)
{
}

CSINSGNSS::CSINSGNSS(int nq0, int nr0, double ts, int yawHkRow0):CSINSTDKF(nq0, nr0)
{
	navStatus = 0;
	posGNSSdelay = vnGNSSdelay = yawGNSSdelay = dtGNSSdelay = dyawGNSS = -0.0f;
	kfts = ts;	gnssLost = &measlost.dd[3];
	lvGNSS = O31;
	Hk(0,3) = Hk(1,4) = Hk(2,5) = 1.0;		// dvn
	Hk(3,6) = Hk(4,7) = Hk(5,8) = 1.0;		// dpos
	yawHkRow = yawHkRow0;
	if(yawHkRow>=6) Hk(yawHkRow,2) = 1.0;	// dyaw
	SetMeasMask(077);
	SetMeasStop(0xffffffff, 1.0);
}

void CSINSGNSS::Init(const CSINS &sins0, int grade)
{
	CSINSTDKF::Init(sins0);
	sins.lever(-lvGNSS, &sins.pos);   // sins0.pos is GNSS pos
	avpi.Init(sins, kfts, 0, 0);
	if(grade==0) {  // inertial grade
		Pmax.Set2(fDEG3(1.0),  fXXX(100.0),  fdPOS(1.0e6), fDPH3(0.5),  fMG3(1.0),
			fXXX(100.0), 0.5, fdKGA15(1000.0,900.0,100.0,100.0));
		Pmin.Set2(fPHI(0.1,1.0),  fXXX(0.001),  fdPOS(.01),  fDPH3(0.001),  fUG3(10.0),
			fXXX(0.001), 0.0001, fdKGA15(1.0,1.0,1.0,1.0));
		Pk.SetDiag2(fPHI(60,600),  fXXX(1.0),  fdPOS(100.0),  fDPH3(0.1),  fMG3(1.0),
			fXXX(1.0),  0.01,  fdKGA15(100.0,90.0,10.0,10.0));
		Qt.Set2(fDPSH3(0.001),  fUGPSHZ3(1.0),  fOOO,  fOO6,
			fOOO, 0.0,  fOO9,  fOO6);
		//MarkovGyro(I31*1000.0, I31*0.01*DPH);  MarkovAcc(I31*1000.0, I31*1.0*UG);
		Xmax.Set(fINF9,  fDPH3(0.5),  fMG3(1.0),
			fXXX(10.0), 0.5,  fdKGA15(1000.0,900.0,1000.0,900.0));
		Rt.Set2(fXXZ(0.5,1.0),   fdLLH(10.0,30.0));  Rt0 = Rt;
		Rmax = Rt*100;  Rmin = Rt*0.01;  Rb = 0.6f;
		innoMax.Set(fXXX(1.0), fdPOS(30.0));  SetInnoMaxcount(5);
		Zmax.Set(fXXX(10.0), fdPOS(1.0e6));
		innoDiffMax.Set(fXXZ(3.0,2.0), fdLLH(100.0,10.0));
		FBOne1.Set(fPHI(1,1), fXXX(0.1), fdLLH(1,3), fDPH3(0.01), fUG3(10), fXXX(0.1), 0.01, fINF9,  fINF6);
		FBTau.Set(fXX9(0.1),  fXX6(1.0),  fINF3, INF,  fINF9,  fINF6);
	}
	else if(grade==1) {  // MEMS grade
		Pmax.Set2(fDEG3(50.0),  fXXX(100.0),  fdPOS(1.0e6), fDPH3(5000.0),  fMG3(30.0),
			fXXX(100.0), 0.5, fdKGA15(1000.0,900.0,100.0,100.0));
		Pmin.Set2(fPHI(1,1),  fXXX(0.0001),  fdPOS(.001),  fDPH3(0.001),  fUG3(1.0),
			fXXX(0.001), 0.0001, fdKGA15(1.0,1.0,1.0,1.0));
		Pk.SetDiag2(fPHI(600,600),  fXXX(1.0),  fdPOS(100.0),  fDPH3(1000.0),  fMG3(10.0),
			fXXX(1.0),  0.01,  fdKGA15(1000.0,90.0,10.0,10.0));
		Qt.Set2(fDPSH3(1.1),  fUGPSHZ3(10.0),  fOOO,  fOO6,
			fOOO, 0.0,  fOO9,  fOO6);
		Xmax.Set(fINF9,  fDPH3(3600.0),  fMG3(50.0),
			fXXX(10.0), 0.5,  fdKGA15(1000.0,900.0,1000.0,900.0));
		Rt.Set2(fXXZ(0.5,1.0),   fdLLH(10.0,30.0));  Rt0 = Rt;
		Rmax = Rt*100;  Rmin = Rt*0.01;  Rb = 0.6f;
		FBTau.Set(fXX9(0.1),  fXX6(1.0),  fINF3, INF,  fINF9,  fINF6);
	}
}

void CSINSGNSS::SetFt(int nnq)
{
	sins.etm();
//	Ft.SetMat3(0,0,sins.Maa), Ft.SetMat3(0,3,sins.Mav), Ft.SetMat3(0,6,sins.Map), Ft.SetMat3(0,9,-sins.Cnb); 
//	Ft.SetMat3(3,0,sins.Mva), Ft.SetMat3(3,3,sins.Mvv), Ft.SetMat3(3,6,sins.Mvp), Ft.SetMat3(3,12,sins.Cnb); 
//						NULL, Ft.SetMat3(6,3,sins.Mpv), Ft.SetMat3(6,6,sins.Mpp);
	Ft.SetMat3(0,0,sins.Maa,sins.Mav,sins.Map), Ft.SetMat3(0,9,-sins.Cnb); 
	Ft.SetMat3(3,0,sins.Mva,sins.Mvv,sins.Mvp), Ft.SetMat3(3,12,sins.Cnb); 
	Ft.SetMat3(6,3,         sins.Mpv,sins.Mpp);
	Ft.SetDiagVect3( 9, 9, sins._betaGyro);
	Ft.SetDiagVect3(12,12, sins._betaAcc);  // 0-14 phi,dvn,dpos,eb,db
	if(nnq==16) {
		Ft(2,15) = -sins.wib.k*sins.Cnb.e22;  // 15 dKGzz
	}
	// if(nnq>=18) NULL;						// 15-17 lever
	// if(nnq>=19) NULL;						// 18 dt
	if(nnq==20) {
		Ft(2,19) = -sins.wib.k*sins.Cnb.e22;  // 19 dKGzz
	}
	else if(nnq==22) {
		CMat3 Cwz=-sins.wib.k*sins.Cnb;
		Ft.SetMat3(0,19, Cwz);				// 19-21 dKG*z
	}
	else if(nnq>=28) {
#ifndef PSINS_FAST_CALCULATION
		CMat3 Cwx=-sins.wib.i*sins.Cnb, Cwy=-sins.wib.j*sins.Cnb, Cwz=-sins.wib.k*sins.Cnb; 
#else
		double wi=-sins.wib.i, wj=-sins.wib.j, wk=-sins.wib.k;
		CMat3 Cwx, Cwy, Cwz; MMULf(Cwx,sins.Cnb,wi); MMULf(Cwy,sins.Cnb,wj); MMULf(Cwz,sins.Cnb,wk);
#endif
//		Ft.SetMat3(0,19, Cwx);  Ft.SetMat3(0,22, Cwy);  Ft.SetMat3(0,25, Cwz);  // 19-27 dKG
		Ft.SetMat3(0,19, Cwx, Cwy, Cwz);  // 19-27 dKG
	}
	if(nnq>=34) {
#ifndef PSINS_FAST_CALCULATION
		CMat3 Cfx= sins.fb.i *sins.Cnb, Cfy= sins.fb.j *sins.Cnb, Cfz= sins.fb.k *sins.Cnb;
#else
		CMat3 Cfx, Cfy, Cfz; MMULf(Cfx,sins.Cnb,sins.fb.i); MMULf(Cfy,sins.Cnb,sins.fb.j); MMULf(Cfz,sins.Cnb,sins.fb.k);
#endif
		Cfz.e00=Cfy.e01,	Cfz.e01=Cfy.e02; 
		Cfz.e10=Cfy.e11,	Cfz.e11=Cfy.e12; 
		Cfz.e20=Cfy.e21,	Cfz.e21=Cfy.e22;
//		Ft.SetMat3(3,28, Cfx);  Ft.SetMat3(3,31, Cfz);  // 28-33 dKA(xx,yx,zx, yy,zy, zz)
		Ft.SetMat3(3,28, Cfx, Cfz);  // 28-33 dKA(xx,yx,zx, yy,zy, zz)
	}
}

void CSINSGNSS::SetHk(int nnq)
{
	if(nnq>=18) {     // GNSS lever
		CMat3 CW=sins.Cnb*askew(sins.webbar), MC=sins.Mpv*sins.Cnb;
		Hk.SetMat3(0,15, -CW);
		Hk.SetMat3(3,15, -MC);
	}
	if(nnq>=19) {    // GNSS dt
		CVect3 MV=sins.Mpv*sins.vn;
//		Hk.SetClmVect3(0,18, -sins.anbar);
		Hk.SetClmVect3(3,18, -MV);
	}
}

void CSINSGNSS::Feedback(int nnq, double fbts)
{
	CKalman::Feedback(nq, fbts);
	for(int i=0; i<avpi.avpinum; i++) {
		avpi.vni[i] -= *(CVect3*)&FBXk.dd[ 3];  avpi.posi[i] -= *(CVect3*)&FBXk.dd[6];
	}
	sins.qnb -= *(CVect3*)&FBXk.dd[0];  sins.vn -= *(CVect3*)&FBXk.dd[ 3];  sins.pos -= *(CVect3*)&FBXk.dd[6];
	sins.eb  += *(CVect3*)&FBXk.dd[9];	sins.db += *(CVect3*)&FBXk.dd[12];  // 0-14 phi,dvn,dpos,eb,db
	// if(nnq==16) {
	// 	double IdKGzz = 1.0-FBXk.dd[15];
	// 	sins.Kg.e20*=IdKGzz, sins.Kg.e21*=IdKGzz, sins.Kg.e22*=IdKGzz;  // 15 dKGzz
	// }
	// if(nnq>=18) {
	// 	lvGNSS += *(CVect3*)&FBXk.dd[15];	// 15-17 lever
	// }
	// if(nnq>=19) {
	// 	dtGNSSdelay += FBXk.dd[18];			// 18 dt
	// }
	// if(nnq==20) {
	// 	double IdKGzz = 1.0-FBXk.dd[19];
	// 	sins.Kg.e20*=IdKGzz, sins.Kg.e21*=IdKGzz, sins.Kg.e22*=IdKGzz;  // 19 dKGzz
	// }
	// else if(nnq==22) {
	// 	CMat3 IdKGz(1.0,0.0,-FBXk.dd[19], 0.0,1.0,-FBXk.dd[20], 0.0,0.0,1.0-FBXk.dd[21]);
	// 	sins.Kg = IdKGz*sins.Kg;		// 19-21 dKG*z
	// }
	// else if(nnq>=28) {
	// 	CMat3 IdKG = I33-(~(*(CMat3*)&FBXk.dd[19]));
	// 	sins.Kg = IdKG*sins.Kg;			// 19-27 dKG
	// }
	// if(nnq>=34) {
	// 	CMat3 IdKA(1.0-FBXk.dd[28],            0.0,           0.0,
	// 		          -FBXk.dd[29],1.0-FBXk.dd[31],           0.0,
	// 				  -FBXk.dd[30],   -FBXk.dd[32],1.0-FBXk.dd[33]);
	// 	sins.Ka = IdKA*sins.Ka;			// 28-33 dKA
	// }
}

void CSINSGNSS::SetMeasGNSS(const CVect3 &posgnss, const CVect3 &vngnss, double yawgnss)
{
	if(!IsZero(posgnss) && avpi.Interp(posGNSSdelay+dtGNSSdelay,0x4))
	{
		*(CVect3*)&Zk.dd[3] = avpi.pos - posgnss;
		SetMeasFlag(00070);
	}
	if(!IsZero(vngnss) && avpi.Interp(vnGNSSdelay+dtGNSSdelay,0x2))
	{
		*(CVect3*)&Zk.dd[0] = avpi.vn - vngnss;
		SetMeasFlag(00007);
	}
	if(!IsZero(yawgnss) && avpi.Interp(yawGNSSdelay+dtGNSSdelay,0x1))
	{
		Zk.dd[yawHkRow] = -diffYaw(avpi.att.k, yawgnss+dyawGNSS);
		SetMeasFlag(01<<yawHkRow);
	}
}

void CSINSGNSS::MeasGNSSZvStop(CVect3 &dvnth, double stop)  // Zv/dvn threshold value
{
	int meas=001;
	double *z=&Zk.dd[0], *v=&dvnth.i;
	for(int i=0; i<3; i++,z++,v++) { 
		if(*z>*v || *z<-(*v))	{ SetMeasStop(meas,stop); meas<<=1; }
	}
}

void CSINSGNSS::MeasGNSSZpStop(CVect3 &dposth, double stop)  // Zk/dpos threshold value
{
	int meas=010;
	double *z=&Zk.dd[3], *p=&dposth.i;
	for(int i=0; i<3; i++,z++,p++) { 
		if(*z>*p || *z<-(*p))	{ SetMeasStop(meas,stop); meas<<=1; }
	}
}

void CSINSGNSS::MeasGNSSZp2X(CVect3 &dposth)  // Zk/dpos threshold value
{
	double *z=&Zk.dd[3], *p=&dposth.i, *ps=&Pset.dd[6];
	for(int i=0; i<3; i++,z++,p++,ps++) {
		if(*z>*p || *z<-(*p))	{ *ps=100.0*(*z)*(*z);	}
	}
}

void CSINSGNSS::Leveling(void)
{
	CVect3 phi=*(CVect3*)&Xk.dd[0]; phi.k=0;
	sins.qnb -= phi;  Xk.dd[0]=Xk.dd[1]=0.0;
	sins.vn -= *(CVect3*)&Xk.dd[3];  *(CVect3*)&Xk.dd[3]=O31;
}

int CSINSGNSS::Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts, int nSteps)
{
	int res=TDUpdate(pwm, pvm, nn, ts, nSteps);
	sins.lever(lvGNSS);
	avpi.Push(sins, 1);
	return res;
}

//***************************  class CEarth  *********************************/
CEarth::CEarth(double a0, double f0, double g0)
{
	a = a0;	f = f0; wie = glv.wie; 
	b = (1-f)*a;
	e = sqrt(a*a-b*b)/a;	e2 = e*e;
	gn = O31;  pgn = 0;  pos0=One31;
	Update(O31);
}

void CEarth::Update(const CVect3 &pos00, const CVect3 &vn0, int isMemsgrade)
{
	this->pos = pos00;  this->vn = vn0;
	CVect3 dpos; VSUB(dpos, pos, pos0);
//	if(dpos.i>10.0/RE||dpos.i<-10.0/RE || dpos.k>1.0||dpos.k<-1.0)  // fast
	{
		sl = sin(pos.i), cl = cos(pos.i), tl = sl/cl;
		double sq = 1-e2*sl*sl, sq2 = sqrt(sq);
		RMh = a*(1-e2)/sq/sq2+pos.k;	f_RMh = 1.0/RMh;
		RNh = a/sq2+pos.k;    clRNh = cl*RNh;  f_RNh = 1.0/RNh; f_clRNh = 1.0/clRNh;
		VEQU(pos0, pos);
	}
	if(isMemsgrade) {
//		wnie.i = 0.0,			wnie.j = wie*cl,		wnie.k = wie*sl;
//		wnen.i = -vn.j*f_RMh,	wnen.j = vn.i*f_RNh,	wnen.k = wnen.j*tl;
		wnin = wnie = wnen = O31;
		sl2 = sl*sl;
		gn.k = -( glv.g0*(1+5.27094e-3*sl2)-3.086e-6*pos.k );
		gcc = pgn ? *pgn : gn;
	}
	else {
		wnie.i = 0.0,			wnie.j = wie*cl,		wnie.k = wie*sl;
		wnen.i = -vn.j*f_RMh,	wnen.j = vn.i*f_RNh,	wnen.k = wnen.j*tl;
//		wnin = wnie + wnen;
		wnin.i=/*wnie.i+*/wnen.i, wnin.j=wnie.j+wnen.j, wnin.k=wnie.k+wnen.k;
		sl2 = sl*sl, sl4 = sl2*sl2;
		gn.k = -( glv.g0*(1+5.27094e-3*sl2+2.32718e-5*sl4)-3.086e-6*pos.k );
		gcc = pgn ? *pgn : gn;
//		gcc -= (wnie+wnin)*vn;
		wnnin.i=wnin.i, wnnin.j=wnie.j+wnin.j, wnnin.k=wnie.k+wnin.k;
		gcc.i -= crosI(wnnin,vn), gcc.j -= crosJ(wnnin,vn), gcc.k -= crosK(wnnin,vn);
	}
}


CVect3 CEarth::vn2dpos(const CVect3 &vn0, double ts) const
{
	return CVect3(vn0.j*f_RMh, vn0.i*f_clRNh, vn0.k)*ts;
}

void CEarth::vn2dpos(CVect3 &dpos, const CVect3 &vn0, double ts) const
{
	dpos.i=vn0.j*f_RMh*ts, dpos.j=vn0.i*f_clRNh*ts, dpos.k=vn0.k*ts;
}

//***************************  class CIMU  *********************************/
CIMU::CIMU(void)
{
	nSamples = 1;
	Reset();
	phim = dvbm = wmm = vmm = swmm = svmm = wm_1 = vm_1 = O31;
	pSf = NULL; pTempArray = NULL; pKga = pgSens = pgSens2 = pgSensX = NULL; pKapn = pKa2 = pTau = NULL; prfu = NULL;
	plv = NULL; pCba = NULL;   tk = tGA = 0.0;
	Sfg = Sfa = One31;  Kg = Ka = Cba = iTCba = I33;
	eb = db = lvx = lvy = lvz = Kapn = Ka2 = Q11 = Q12 = Q13 = Q21 = Q22 = Q23 = Q31 = Q32 = Q33 = O31;  SSx = SSy = SSz = O33;
}

void CIMU::Reset(void)
{
	preFirst = onePlusPre = preWb = true;  swmm = svmm = wm_1 = vm_1 = O31;  smmT = 0.0;
	compensated = false;
}

void CIMU::SetRFU(const char *rfu0)
{
	for(int i=0; i<3; i++) rfu[i]=rfu0[i];
	prfu = rfu;
}

void CIMU::SetSf(const CVect3 &Sfg0, const CVect3 &Sfa0)
{
	Sfg = Sfg0;  Sfa = Sfa0;
	pSf = &Sfg;
}

void CIMU::SetTemp(double *tempArray0, int type)
{
	pTempArray = tempArray0;  iTemp = 0;
	double *p=&Kg.e00, *p1=tempArray0;
	for(int k=0; k<37; k++,p++,p1+=5) *p = *p1;
	pKga = &Kg;  if(type>1) pKa2 = &Ka2;  if(type>2) plv = &lvx;
}

void CIMU::SetKga(const CMat3 &Kg0, const CVect3 eb0, const CMat3 &Ka0, const CVect3 &db0)
{
	Kg = Kg0; eb = eb0; Ka = Ka0; db = db0;
	pKga = &Kg;
}

void CIMU::SetgSens(const CMat3 &gSens0, const CMat3 &gSens20, const CMat3 &gSensX0)
{
	gSens = gSens0;	pgSens = &gSens;
	if(&gSens20!=&O33) { gSens2 = gSens20;	pgSens2 = &gSens2; }
	if(&gSensX0!=&O33) { gSensX = gSensX0;	pgSensX = &gSensX; }
}

void CIMU::SetKa2(const CVect3 &Ka20)
{
	Ka2 = Ka20;	pKa2 = &Ka2;
}

void CIMU::SetKapn(const CVect3 &Kapn0)
{
	Kapn = Kapn0;	pKapn = &Kapn;
}

void CIMU::SetLvtGA(const CVect3 &lvx0, const CVect3 &lvy0, const CVect3 &lvz0, double tGA0)
{
	lvx = lvx0; lvy = lvy0; lvz = lvz0; tGA = tGA0;
	plv = &lvx;
}

void CIMU::SetTau(const CVect3 &taux0, const CVect3 &tauy0, const CVect3 &tauz0)
{
	Taux = taux0;	Tauy = tauy0;	Tauz = tauz0;	pTau = &Taux;
}

void CIMU::SetCba(const CMat3 &Cba0)
{
	Cba = Cba0;  iTCba = inv(~Cba);  pCba = &Cba;
	CMat3 U = inv(~Cba);
	CVect3 V1=Cba.GetClm(0), V2=Cba.GetClm(1), V3=Cba.GetClm(2);
    Q11 = U.e00*V1, Q12 = U.e01*V2, Q13 = U.e02*V3,
    Q21 = U.e10*V1, Q22 = U.e11*V2, Q23 = U.e12*V3,
    Q31 = U.e20*V1, Q32 = U.e21*V2, Q33 = U.e22*V3;
}

double CIMU::GetMeanwf(CVect3 &wib, CVect3 &fsf, BOOL reset)
{
	double mT=smmT;
	wib = swmm*(1.0/smmT);  fsf = svmm*(1.0/smmT);
	if(reset) { swmm=svmm=O31; smmT=0.0; }
	return mT;
}

double CIMU::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples0, double ts)
{
	static double conefactors[5][4] = {				// coning coefficients
		{2./3},										// 2
		{9./20, 27./20},							// 3
		{54./105, 92./105, 214./105},				// 4
		{250./504, 525./504, 650./504, 1375./504}	// 5
		};
	int i;
	double *pcf = conefactors[nSamples0-2];
	CVect3 cm(0.0), sm(0.0), wm[6], vm[6], vmma(0.0); wmm=O31, vmm=O31;

	psinsassert(nSamples0>0 && nSamples0<6);
	for(i=0; i<nSamples0; i++) { wm[i]=pwm[i]; vm[i]=pvm[i]; }
	if(pSf)
	{
		for(i=0; i<nSamples0; i++) { wm[i] = dotdiv(wm[i], Sfg);  vm[i] = dotdiv(vm[i], Sfa); }
	}
	if(pKa2||pKapn)
	{
		for(i=0; i<nSamples0; i++) { vmma+=vm[i]; }
	}
	if(pKga)
	{
		for(i=0; i<nSamples0; i++) { wm[i] = Kg*wm[i] - eb*ts;  vm[i] = Ka*vm[i] - db*ts;  }
	}
	this->nSamples = nSamples0; nts=nSamples*ts; _nts=1.0/nts; tk+=nts;
	if(nSamples==1 && onePlusPre)  // one-plus-previous sample
	{
		if(preFirst) { wm_1=wm[0]; vm_1=vm[0]; preFirst=false; }
		cm = 1.0/12*wm_1;
		sm = 1.0/12*vm_1;
	}
	for(i=0; i<nSamples-1; i++)
	{
		cm += pcf[i]*wm[i];
		sm += pcf[i]*vm[i];
		wmm += wm[i];
		vmm += vm[i];
	}
	wm_1=wm[i];  vm_1=vm[i];
	wmm += wm[i];
	vmm += vm[i];
	phim = wmm + cm*wm[i];
	dvbm = vmm + 1.0/2*wmm*vmm + (cm*vm[i]+sm*wm[i]);
	if(nSamples==1 && compensated) 
	{
		phim = wmm;
		dvbm = vmm;
	}
	CVect3 _phim(0.0), _dvbm(0.0);
	if(pTempArray)
	{
		double *p=&pTempArray[iTemp*5];
		*(&Kg.e00+iTemp) = p[0] + polyval(&p[1], 3, Temp);
		if(++iTemp==40) iTemp=0;  // (&tGA-&Kg.e00)==37
	}
	if(pgSens) {
		_phim.i += gSens.e00*dvbm.i+gSens.e01*dvbm.j+gSens.e02*dvbm.k;   // gSens.eij in (rad/s)/(m/ss)
		_phim.j += gSens.e10*dvbm.i+gSens.e11*dvbm.j+gSens.e12*dvbm.k;
		_phim.k += gSens.e20*dvbm.i+gSens.e21*dvbm.j+gSens.e22*dvbm.k;
	}
	if(pgSens2) {
		double fx2_Ts = dvbm.i*dvbm.i*_nts,  fy2_Ts = dvbm.j*dvbm.j*_nts,  fz2_Ts = dvbm.k*dvbm.k*_nts;
		_phim.i += gSens2.e00*fx2_Ts+gSens2.e01*fy2_Ts+gSens2.e02*fz2_Ts;   // gSens2.eij in (rad/s)/(m/ss)^2
		_phim.j += gSens2.e10*fx2_Ts+gSens2.e11*fy2_Ts+gSens2.e12*fz2_Ts;
		_phim.k += gSens2.e20*fx2_Ts+gSens2.e21*fy2_Ts+gSens2.e22*fz2_Ts;
	}
	if(pgSensX) {
		double fxy_Ts = dvbm.i*dvbm.j*_nts,  fyz_Ts = dvbm.j*dvbm.k*_nts,  fzx_Ts = dvbm.k*dvbm.i*_nts;
		_phim.i += gSensX.e00*fxy_Ts+gSensX.e01*fyz_Ts+gSensX.e02*fzx_Ts;   // gSensX.eij in (rad/s)/(m/ss)^2
		_phim.j += gSensX.e10*fxy_Ts+gSensX.e11*fyz_Ts+gSensX.e12*fzx_Ts;
		_phim.k += gSensX.e20*fxy_Ts+gSensX.e21*fyz_Ts+gSensX.e22*fzx_Ts;
	}
	if(pTau) {
//		_phim.i += (dvbm.i*Taux.i+dvbm.j*Taux.j+dvbm.k*Taux.k)*phim.i*_nts;
//		_phim.j += (dvbm.i*Tauy.i+dvbm.j*Tauy.j+dvbm.k*Tauy.k)*phim.j*_nts;
//		_phim.k += (dvbm.i*Tauz.i+dvbm.j*Tauz.j+dvbm.k*Tauz.k)*phim.k*_nts;
		CVect3 wXf=wmm*vmm;
		_phim.i += dot(Taux, wXf)*_nts;
		_phim.j += dot(Tauy, wXf)*_nts;
		_phim.k += dot(Tauz, wXf)*_nts;
	}
	if(pKapn) {
		_dvbm += iTCba*dotmul(Kapn,abs(vmma));  // Kapn in (m/ss)/(m/ss)=(m/s)/(m/s)=no unit
	}
	if(pKa2) {
		vmma.i = Ka2.i*vmma.i*vmma.i*_nts;   // Ka2.i in (m/ss)/(m/ss)^2
		vmma.j = Ka2.j*vmma.j*vmma.j*_nts;
		vmma.k = Ka2.k*vmma.k*vmma.k*_nts;
		_dvbm += iTCba*vmma;
	}
	if(plv) {
		if(preWb) { wb_1=phim*_nts, preWb=false; }
		CVect3 wb=phim*_nts, dwb=(wb-wb_1)*_nts;  wb_1=wb;
		CMat3 W=askew(dwb)+pow(askew(wb),2);
		CVect3 fL;
		if(pCba) {
			SSx = CMat3(Q11*W, Q21*W, Q31*W);  SSy = CMat3(Q12*W, Q22*W, Q32*W);  SSz = CMat3(Q13*W, Q23*W, Q33*W);
			fL = SSx*lvx + SSy*lvy + SSz*lvz;
		}
		else {
			fL = CVect3(dot(*(CVect3*)&W.e00,lvx), dot(*(CVect3*)&W.e10,lvy), dot(*(CVect3*)&W.e20,lvz));
		}
		_dvbm += fL*nts + tGA*(wb*dvbm);
	}
	wmm -= _phim;	  vmm -= _dvbm;
	phim -= _phim;	  dvbm -= _dvbm;
	swmm += wmm;     svmm += vmm;  smmT += nts;
	if(prfu) { IMURFU(&wmm, &vmm, 1, prfu);  IMURFU(&phim, &dvbm, 1, prfu); IMURFU(&swmm, &svmm, 1, prfu); }
	return tk;
}

void IMURFU(CVect3 *pwm, int nSamples, const char *str)
{
	if(str[0]=='X') return;
	for(int n=0; n<nSamples; n++)
	{
		CVect3 tmpwm;
		double *pw=(double*)&pwm[n].i;
		for(int i=0; i<3; i++,pw++)
		{
			switch(str[i])
			{
			case 'R':  tmpwm.i= *pw;  break;
			case 'L':  tmpwm.i=-*pw;  break;
			case 'F':  tmpwm.j= *pw;  break;
			case 'B':  tmpwm.j=-*pw;  break;
			case 'U':  tmpwm.k= *pw;  break;
			case 'D':  tmpwm.k=-*pw;  break;
			}
		}
		pwm[n] = tmpwm;
	}
}

void IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str)
{
	if(str[0]=='X') return;
	IMURFU(pwm, nSamples, str);
	IMURFU(pvm, nSamples, str);
}

void IMUStatic(CVect3 &wm, CVect3 &vm, CVect3 &att0, CVect3 &pos0, double ts)
{
	CEarth eth;
	eth.Update(pos0);
	CMat3 Cbn=~a2mat(att0);
	wm = Cbn*eth.wnie*ts;	vm = -Cbn*eth.gn*ts;
}

CMat lsclbt(CMat &wfb, CMat &wfn)
{
	CMat A(wfb.row, 4, 0.0), Kga_edb(4,3,0.0);
	CVect Y(wfb.row), X;
	for(int k=0; k<3; k++)
	{
		for(int i=0; i<wfb.row; i++)
		{
			A.SetRow(i, CVect(4, wfb(i,0), wfb(i,1), wfb(i,2), -1.0));
			Y(i) = wfn(i,k);
		}
		X = inv4((~A)*A)*((~A)*Y);
		Kga_edb(k,0)=X(0), Kga_edb(k,1)=X(1), Kga_edb(k,2)=X(2), Kga_edb(3,k)=X(3);
	}
	return Kga_edb;
}

//**************************  class CIMUInc  ********************************/
CIMUInc::CIMUInc(double gScale0, double aScale0)
{
	Init();
}

void CIMUInc::Init(double gScale0, double aScale0)
{
	diGx=diGy=diGz=diAx=diAy=diAz=0;
	iTotalGx0=iTotlaGy0=iTotalGz0=iTotalAx0=iTotalAy0=iTotalAz0=0;
	iTotalGx=iTotalGy=iTotalGz=iTotalAx=iTotalAy=iTotalAz=0;
	fTotalGx=fTotalGy=fTotalGz=fTotalAx=fTotalAy=fTotalAz=0.0;
	this->gScale = gScale0,  this->aScale = aScale0;
}

void CIMUInc::Update(const CVect3 &wm, const CVect3 &vm)
{
#define MAX_NN  2147483648.0  // 2^31
#define MAX_PP  (MAX_NN-1.0)
	fTotalGx += wm.i/gScale, fTotalGy += wm.j/gScale, fTotalGz += wm.k/gScale;
	fTotalAx += vm.i/aScale, fTotalAy += vm.j/aScale, fTotalAz += vm.k/aScale;
	double *pfTotal=&fTotalGx;  int *pdi = &diGx, *piTotal0 = &iTotalGx0, *piTotal = &iTotalGx;
	for(int i=0; i<6; i++,pfTotal++,pdi++,piTotal0++,piTotal++) {
		if(*pfTotal>=MAX_PP)		*pfTotal -= 2.0*MAX_NN;
		else if(*pfTotal<=-MAX_NN) 	*pfTotal += 2.0*MAX_NN;
		*piTotal = (int)(*pfTotal);
		*pdi = *piTotal - *piTotal0;
		if(*pdi>100000000)			*pdi -= 2147483648;
		else if(*pdi<-100000000)	*pdi += 2147483648;
		*piTotal0 = *piTotal;
	}
}

//**************************  class CIMUInv  ********************************/
CIMUInv::CIMUInv(const CVect3 &att00, const CVect3 &vn00, const CVect3 &pos00, double ts0, double tk0)
{
	vn0 = vn00;  pos0 = pos00;  att = att00;  vn = vn0;  pos = pos0;
	Cbn0 = ~a2mat(att00);
	eth.Update(pos, vn);
	ts = ts0;  tk = tk0;
	isFirst = TRUE;
}

void CIMUInv::Update(const CVect3 &att1, const CVect3 &pos1)
{
	CVect3 vn1 = pp2vn(pos1, pos0, ts, &eth);  vn1 = vn0+(vn1-vn0)*2;
	eth.Update((pos0+pos1)/2, (vn0+vn)/2);
	CMat3 Cnb1 = a2mat(att1);
	CVect3 phim = m2rv(Cbn0*rv2m(eth.wnin*ts)*Cnb1);
	if(isFirst) wm0=phim;
    wm = inv(I33+askew(wm0/12.0))*phim;
    CVect3 dvbm = Cbn0 * (rv2q(eth.wnin*ts*2) * (vn1-vn0-eth.gcc*ts));
	if(isFirst) { vm0=dvbm; isFirst=FALSE; }
    vm = inv(I33+askew(wm/2.0+wm0/12.0))*(dvbm-vm0*wm/12.0);
	wm0 = wm;  vm0 = vm;
	Cbn0 = ~Cnb1;  vn0 = vn1;  pos0 = pos1;
	att = att1;  vn = vn1;  pos = pos1;
	tk += ts;
}

//***************************  class CSINS  *********************************/
CSINS::CSINS(double yaw0, const CVect3 &pos0, double tk0)
{
	Init(a2qua(CVect3(0,0,yaw0)), O31, pos0, tk0);
}

CSINS::CSINS(const CVect3 &att0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
	Init(a2qua(att0), vn0, pos0, tk0);
}

CSINS::CSINS(const CQuat &qnb0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
	Init(qnb0, vn0, pos0, tk0);
}

void CSINS::Init(const CVect3 &att0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
	Init(a2qua(att0), vn0, pos0, tk0);
}

void CSINS::Init(const CQuat &qnb0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
	tk = tk0;  ts = nts = 1.0;  dist = -EPS;
	velMax = 400.0; hgtMin = -RE*0.01, hgtMax = -hgtMin; latMax = 85.0*DEG; afabar = 0.1;
	qnb = qnb0;	vn = vn0, pos = pos0;
	Kg = Ka = I33; eb = db = Ka2 = O31;
	Maa = Mav = Map = Mva = Mvv = Mvp = Mpv = Mpp = O33;
	SetTauGA(CVect3(INF),CVect3(INF));
	CVect3 wib0(0.0), fb0=(~qnb)*CVect3(0,0,glv.g0);
	lvr = an = anbar = webbar = O31;
	isOpenloop = isMemsgrade = isNocompasseffect = isOutlever = 0;
	Cbn = I33;
	Update(&wib0, &fb0, 1, 1.0); imu.preFirst = 1;
	tk = tk0;  ts = nts = 1.0; qnb = qnb0;	att=q2att(qnb), vn = vn0, pos = pos0;
	mmwb = mmfb = CMaxMin(100);
	mvn = mvnmax = mvnmin = vn; mvni = O31; mvnt = mvnT = lvlT = 0.0; mvnk = 0;
	etm(); lever(); Extrap();
}

void CSINS::SetTauGA(const CVect3 &tauG, const CVect3 &tauA)
{
	if(tauG.i>EPS) {
		tauGyro = tauG;
		_betaGyro.i = tauG.i>INFp5 ? 0.0 : -1.0/tauG.i;   // Gyro&Acc inverse correlation time for AR(1) model
		_betaGyro.j = tauG.j>INFp5 ? 0.0 : -1.0/tauG.j;
		_betaGyro.k = tauG.k>INFp5 ? 0.0 : -1.0/tauG.k;
	}
	if(tauA.i>EPS) {
		tauAcc = tauA;
		_betaAcc.i  = tauA.i>INFp5 ? 0.0 : -1.0/tauA.i;
		_betaAcc.j  = tauA.j>INFp5 ? 0.0 : -1.0/tauA.j;
		_betaAcc.k  = tauA.k>INFp5 ? 0.0 : -1.0/tauA.k;
	}
}

void CSINS::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts0)
{
	this->ts = ts0;  nts = nSamples*ts;	tk += nts;
	double nts2 = nts/2, _nts=1.0/nts;
	if(isMemsgrade) {
		imu.Update(pwm, pvm, nSamples, ts);
		imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
		if(!isOpenloop) eth.Update(pos,O31,1);
		wib = imu.phim*_nts; fb = imu.dvbm*_nts;
		web = wib;  webbar = (1-afabar)*webbar + afabar*web;
		wnb = wib;
		fn = qnb*fb;
		an = fn+eth.gcc;  anbar = (1-afabar)*anbar + afabar*an;
		CVect3 vn1 = vn + an*nts;
		pos = pos + eth.vn2dpos(vn+vn1, nts2);	vn = vn1;
		qnb = qnb*rv2q(imu.phim);
		Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
	}
	else {
#ifndef PSINS_FAST_CALCULATION
		imu.Update(pwm, pvm, nSamples, ts);
		imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
		CVect3 vn01 = vn+an*nts2, pos01 = pos+eth.vn2dpos(vn01,nts2);
		if(!isOpenloop) eth.Update(pos01, vn01);
		wib = imu.phim/nts; fb = imu.dvbm/nts;
		web = wib - Cbn*eth.wnie;  webbar = (1-afabar)*webbar + afabar*web;
		wnb = wib - Cbn*eth.wnin;
		fn = qnb*fb;
		an = rv2q(-eth.wnin*nts2)*fn+eth.gcc;  anbar = (1-afabar)*anbar + afabar*an;
		CVect3 vn1 = vn + an*nts;
		pos = pos + eth.vn2dpos(vn+vn1, nts2);	vn = vn1;
		qnb = rv2q(-eth.wnin*nts)*qnb*rv2q(imu.phim);
		Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
#else  // fast
		imu.Update(pwm, pvm, nSamples, ts);
//		imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
		AXbt(imu.phim, Kg, imu.phim, eb, -nts);  AXbt(imu.dvbm, Ka, imu.dvbm, db, -nts);
//		CVect3 vn01 = vn+an*nts2, pos01 = pos+eth.vn2dpos(vn01,nts2);
		CVect3 vn01,pos01,dpos; VADDf(vn01,vn,an,nts2); eth.vn2dpos(dpos,vn01,nts2); VADD(pos01,pos,dpos);
		if(!isOpenloop) eth.Update(pos01, vn01);
//		wib = imu.phim*_nts; fb = imu.dvbm*_nts;
		VMULf(wib, imu.phim, _nts);  VMULf(fb, imu.dvbm, _nts);
//		web = wib - Cbn*eth.wnie;  webbar = (1-afabar)*webbar + afabar*web;
		CVect3 vtmp; mul(vtmp, Cbn, eth.wnie); VSUB(web, wib, vtmp); double b=1-afabar; VADDff(webbar, webbar, b, web, afabar);
//		wnb = wib - Cbn*eth.wnin;
		mul(vtmp, Cbn, eth.wnin); VSUB(wnb, wib, vtmp);
//		fn = qnb*fb;
		mul(fn, qnb, fb);
//		an = rv2q(-eth.wnin*nts2)*fn+eth.gcc;  anbar = (1-afabar)*anbar + afabar*an;
		CVect3 wnint, fnn; double tn=-nts2; VMULf(wnint,eth.wnin,tn);
		CQuat q; rv2q(q, wnint); mul(fnn, q, fn);
		VADD(an, fnn, eth.gcc);  VADDff(anbar, anbar, b, an, afabar);
//		CVect3 vn1 = vn + an*nts;
		CVect3 vn1; VADDf(vn1, vn, an, nts);
//		pos = pos + eth.vn2dpos(vn+vn1, nts2);	vn = vn1;
		VADD(vn01, vn, vn1);  eth.vn2dpos(dpos, vn01, nts2);  VADD(pos, pos, dpos);  VEQU(vn, vn1);
//		qnb = rv2q(-eth.wnin*nts)*qnb*rv2q(imu.phim);
		tn=-nts; VMULf(wnint,eth.wnin,tn);	rv2q(q, wnint);  mul(qnb, q, qnb);  rv2q(q, imu.phim);  mul(qnb, qnb, q); 
//		Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
		q2mat(Cnb, qnb);  m2att(att, Cnb);  _TT(Cbn, Cnb);  mul(vb, Cbn, vn);
#endif		
	}
	psinsassert(pos.i<85.0*DEG && pos.i>-85.0*DEG);
	if(vn.i>velMax) vn.i=velMax; else if(vn.i<-velMax) vn.i=-velMax;
	if(vn.j>velMax) vn.j=velMax; else if(vn.j<-velMax) vn.j=-velMax;
	if(vn.k>velMax) vn.k=velMax; else if(vn.k<-velMax) vn.k=-velMax;
	if(pos.i>latMax) pos.i=latMax; else if(pos.i<-latMax) pos.i=-latMax;
	if(pos.j>PI) pos.j-=_2PI; else if(pos.j<-PI) pos.j+=_2PI;
	if(pos.k>hgtMax) pos.k=hgtMax; else if(pos.k<hgtMin) pos.k=hgtMin;
	if(mvnT>EPS) {   // calculate max/min/mean-vn within mvnT
		if(mvnk==0) {
			mvnmax = mvnmin = vn;
			mvni = O31;  mvnt = 0.0;  mvnCnb0 = Cnb;
		}
		else {
			if(vn.i>mvnmax.i) mvnmax.i=vn.i; else if(vn.i<mvnmin.i) mvnmin.i=vn.i; 
			if(vn.j>mvnmax.j) mvnmax.j=vn.j; else if(vn.j<mvnmin.j) mvnmin.j=vn.j; 
			if(vn.k>mvnmax.k) mvnmax.k=vn.k; else if(vn.k<mvnmin.k) mvnmin.k=vn.k; 
		}
		mvni += vn;  mvnt += nts;  mvnk++;
		if(mvnt>=mvnT) {
			mvn = mvni*(1.0/mvnk);
			mvnk = 0;   // OK if mvnk==0
		}
	}
	if(dist>=0.0)
		dist += sqrt(vn.i*vn.i + vn.j*vn.j + vn.k*vn.k)*nts;  // travel distance ^2
	if(lvlT>=0)
		lvlT += nts;
	mmwb.Update((float)(norm(imu.wmm)*_nts));
	mmfb.Update((float)(norm(imu.vmm)*_nts+eth.gn.k));
}

void CSINS::Extrap(const CVect3 &wm, const CVect3 &vm, double ts0)
{
	if(ts0<1.0e-6)  // reset
	{
		qnbE = qnb, vnE = vn, posE = pos, attE = att;
	}
	else
	{
		vnE = vnE + qnbE*vm + eth.gcc*ts0;
		posE = posE + eth.vn2dpos(vnE,ts0);
		qnbE = qnbE*rv2q(wm); attE = q2att(qnbE);
	}
}

void CSINS::Extrap(double extts)
{
	double k = extts/nts;
	vnE = vn + qnb*imu.dvbm*k + eth.gcc*extts;
	posE = pos + eth.vn2dpos(vn,extts);
	attE = q2att(qnb*rv2q(imu.phim*k));
}

void CSINS::lever(const CVect3 &dL, CVect3 *ppos, CVect3 *pvn)
{
	if(&dL!=&O31) lvr = dL;
//	Mpv = CMat3(0,eth.f_RMh,0, eth.f_clRNh,0,0, 0,0,1);
	Mpv.e01=eth.f_RMh, Mpv.e10=eth.f_clRNh, Mpv.e22=1.0;
	CW = Cnb*askew(web), MpvCnb = Mpv*Cnb;
	if(ppos==NULL) {
		posL = pos + MpvCnb*lvr;
		if(pvn==NULL)  vnL = vn + CW*lvr;
	} 
	else {
		*ppos = pos + MpvCnb*lvr;
		if(pvn!=NULL)  *pvn = vn + CW*lvr;
	}
}

void CSINS::lever2(const CVect3 &dL, CVect3 *ppos, CVect3 *pvn, const CVect3 *ppos0, const CVect3 *pvn0)
{
	if(&dL!=&O31) lvr = dL;
	if(ppos0==NULL) ppos0=ppos;
	if(pvn0==NULL) pvn0=pvn;
	Mpv.e01=eth.f_RMh, Mpv.e10=eth.f_clRNh, Mpv.e22=1.0;
	CW = Cnb*askew(web), MpvCnb = Mpv*Cnb;
	if(ppos!=NULL)  *ppos = *ppos0 + MpvCnb*lvr;
	if(pvn!=NULL)  *pvn = *pvn0 + CW*lvr;
}

void CSINS::atss(double *attack, double *sideslip)
{
	CVect3 as=::atss(att, vn);
	*attack = as.i,  *sideslip = as.k;
}

void CSINS::etm(void)
{
	if(isMemsgrade) {
		Mva = askew(fn);
		Mpv = CMat3(0,eth.f_RMh,0, eth.f_clRNh,0,0, 0,0,1);
	}
	else {
		double tl=eth.tl, secl=1.0/eth.cl, secl2=secl*secl, 
			wN=eth.wnie.j, wU=eth.wnie.k, vE=vn.i, vN=vn.j;
		double f_RMh=eth.f_RMh, f_RNh=eth.f_RNh, f_clRNh=eth.f_clRNh, 
			f_RMh2=f_RMh*f_RMh, f_RNh2=f_RNh*f_RNh;
		CMat3 Avn=askew(vn),
			Mp1(0,0,0, -wU,0,0, wN,0,0),
			Mp2(0,0,vN*f_RMh2, 0,0,-vE*f_RNh2, vE*secl2*f_RNh,0,-vE*tl*f_RNh2);
		if(isNocompasseffect)	Maa = O33;
		else					Maa = askew(-eth.wnin);
		Mav = CMat3(0,-f_RMh,0, f_RNh,0,0, tl*f_RNh,0,0);
		Map = Mp1+Mp2;
		Mva = askew(fn);
		Mvv = Avn*Mav - askew(eth.wnie+eth.wnin);
		Mvp = Avn*(Mp1+Map);
		double scl = eth.sl*eth.cl;
		Mvp.e20 = Mvp.e20-glv.g0*(5.27094e-3*2*scl+2.32718e-5*4*eth.sl2*scl); Mvp.e22 = Mvp.e22+3.086e-6;
		Mpv = CMat3(0,f_RMh,0, f_clRNh,0,0, 0,0,1);
		Mpp = CMat3(0,0,-vN*f_RMh2, vE*tl*f_clRNh,0,-vE*secl*f_RNh2, 0,0,0);
	}
}

void CSINS::AddErr(const CVect3 &phi, const CVect3 &dvn, const CVect3 &dpos)
{
	qnb -= -phi;
	vn += dvn;
	pos += CVect3(dpos.i*eth.f_RMh, dpos.j*eth.f_clRNh, dpos.k);  // NOTE: dpos in meter
}

void CSINS::AddErr(double phiU, const CVect3 &dvn, const CVect3 &dpos)
{
	AddErr(CVect3(0,0,phiU), dvn, dpos);
}

void CSINS::Leveling(int flag)
{
	if(flag==0) {
		lvlVn0 = vn;  lvlT = EPS;  // record
	}
	else if(flag>=1)
	{
		CVect3 dvn = vn-lvlVn0;  vn = O31;
		CVect3 phi = CVect3(dvn.j, -dvn.i, 0)/(lvlT*G0);
		qnb -= phi;   // rectify
		if(flag==2)  db.k += dvn.k/lvlT;
		lvlT = -1.0;  // stop
	}
}

void CSINS::DebugStop(double t1, int absT, int ext)
{
	static double dst0=-INF;
	if(dst0<-INFp5) dst0=tk;
	double t=tk-dst0*absT;  // absT=0 for absolute time
	if(t1<t && t<t1+3*nts) {
		if(ext==0)
			t1 = tk-dst0*absT;  // please place breakpoint here
		else
			exit(0);  // for exit
	}
}

//***************************  class CDR  *********************************/
CDR::CDR(void)
{
}

void CDR::Init(const CSINS &sins, const CVect3 &kappa)
{
	Init(sins.att, sins.pos, kappa, sins.tk);
	vn = sins.vn;  velPre = norm(vn);
}

void CDR::Init(const CVect3 &att0, const CVect3 &pos0, const CVect3 &kappa, double tk0)
{
	Kod = kappa.j;
	Cbo = a2mat(CVect3(kappa.i,0.0,kappa.k));
	qnb = a2qua(att0);  Cnb = q2mat(qnb); att = m2att(Cnb);
	vn = O31;	pos = pos0;  eth.Update(pos, vn);
	velPre = 0.0;  velMax=50.0; velMin=-3.0;  afa = 0.1;
	tk = tk0;  dist = 0.0;
	SetGCK(10.0);
}

void CDR::Update(const CVect3 &wm, double dS, double ts, const CVect3 &vm)
{
	tk += ts;  dist += dS;
	double vel = dS*Kod/ts;
	if(vel>velMax||vel<velMin) vel=velPre;  // if abnormal, keep last velocity
	CVect3 vnk = Cnb*(Cbo*CVect3(0,vel,0));  velPre=vel;
	eth.Update(pos, vnk);
	vn = (1-afa)*vn + afa*vnk;  // AR(1) vel filter
	pos = pos + eth.vn2dpos(vnk, ts); 
	qnb = rv2q(-eth.wnin*ts)*qnb*rv2q(wm);
	Cnb = q2mat(qnb); att = m2att(Cnb);
	if(&vm!=&O31) {
		Leveling(vm, ts);
	}
}

void CDR::SetGCK(double Td)
{
	double xi=0.707, xi2=xi*xi, ws2=G0/RE, sigma=_2PI*xi/(Td*sqrt(1.0-xi2)), sigma2=sigma*sigma;
    gck1 = 3.0*sigma; 
    gck2 = sigma2*(2.0+1.0/xi2)/ws2-1.0; 
    gck3 = sigma2*sigma/(G0*xi2);
	wnc = vni = dpos = O31;
}

void CDR::Leveling(const CVect3 &vm, double ts)
{
	CVect3 fn=qnb*vm/ts;
	double dVE = vni.i - vn.i;     // vni: inertial vel;  vn: ref vel
	vni.i = vni.i + (fn.i-gck1*dVE)*ts;
	dpos.i = dpos.i + dVE*gck3*ts;
	wnc.j = dVE*(1+gck2)/RE + dpos.i;
	double dVN = vni.j - vn.j;
	vni.j = vni.j + (fn.j-gck1*dVN)*ts;
	dpos.j = dpos.j + dVN*gck3*ts;
	wnc.i = -dVN*(1+gck2)/RE - dpos.j;
	qnb = rv2q(-wnc*ts)*qnb;
}

CVect3 CDR::Calibrate(const CVect3 &pos0, const CVect3 &pos1, const CVect3 &pos1DR, double dist0)
{
	CVect3 dpos0=pos1-pos0, xyz, xyzDR;
	xyz = CVect3(dpos0.j*eth.clRNh, dpos0.i*eth.RMh, dpos0.k);
	dpos0=pos1DR-pos0;
	xyzDR = CVect3(dpos0.j*eth.clRNh, dpos0.i*eth.RMh, dpos0.k);
	if(dist0<-EPS) dist0=-norm(xyz);
	else if(dist0<1.0) dist0=norm(xyz);
	double kod = norm(xyz)/norm(xyzDR);
	double dpitch = (xyzDR.k-xyz.k)/dist0;
	double dyaw=crossXY(xyz,xyzDR)/normXY(xyz)/normXY(xyzDR);
	return CVect3(dpitch, kod, dyaw);
}

//***************************  class CAVPInterp  *********************************/
CAVPInterp::CAVPInterp(void)
{
}

void CAVPInterp::Init(const CSINS &sins, double ts0, BOOL islever, int num)
{
	if(islever)
		Init(sins.att, sins.vnL, sins.posL, ts0, num);
	else
		Init(sins.att, sins.vn, sins.pos, ts0, num);
}

void CAVPInterp::Init(const CVect3 &att0, const CVect3 &vn0, const CVect3 &pos0, double ts0, int num)
{
	psinsassert(num<AVPINUM);
	this->ts = ts0;
	ipush = 0;  avpinum = num;
	for(int i=0; i<AVPINUM; i++) { atti[i]=att0, vni[i]=vn0; posi[i]=pos0; }
	att = att0, vn = vn0, pos = pos0;
}

void CAVPInterp::Push(const CSINS &sins, BOOL islever)
{
	if(islever)
		Push(sins.att, sins.vnL, sins.posL);
	else
		Push(sins.att, sins.vn, sins.pos);
}

void CAVPInterp::Push(const CVect3 &attk, const CVect3 &vnk, const CVect3 &posk)
{
	if(++ipush>=avpinum) ipush = 0;
	atti[ipush] = attk; vni[ipush] = vnk; posi[ipush] = posk;
}

int CAVPInterp::Interp(double tpast, int avp)
{
	int res=1, k, k1, k2;
	if(tpast<-avpinum*ts) tpast=-avpinum*ts; else if(tpast>0) tpast=0;
//	if(tpast<-AVPINUM*ts||tpast>0) return (res=0);
	k = (int)(-tpast/ts);
	if((k2=ipush-k)<0) k2 += avpinum;
	if((k1=k2-1)<0)  k1 += avpinum;
	double tleft = -tpast - k*ts;
	if(tleft>0.99*ts)	{
		if(avp&0x1) att=atti[k1];
		if(avp&0x2) vn=vni[k1]; 
		if(avp&0x4) pos=posi[k1];
	}
	else if(tleft<0.01*ts)	{
		if(avp&0x1) att=atti[k2]; 
		if(avp&0x2) vn=vni[k2]; 
		if(avp&0x4) pos=posi[k2];
	}
	else	{
		double b=tleft/ts, a=1-b;
		if(avp&0x1)	
		{ 
			att = b*atti[k1]+a*atti[k2]; 
			if(normInf(att-atti[k1])>10.0*DEG) res=0;
		}
		if(avp&0x2)	{ vn  = b*vni[k1] +a*vni[k2]; }
		if(avp&0x4)	
		{
			pos = b*posi[k1]+a*posi[k2];
			if(fabs(pos.j-posi[k1].j)>1.0*DEG) res=0;
		}
	}
	return res;
}


BYTE* flipud(BYTE *p, int rows, int clmBytes)
{
	BYTE *ptmp=(BYTE*)malloc(clmBytes), *p0=p, *p1=p0+(rows-1)*clmBytes, *p11=p1;
	memcpy(ptmp, p0, clmBytes);
	while(1) {
		memcpy(p0, p1, clmBytes);  p0+=clmBytes;
		if(p0>=p1) break;
		memcpy(p1, p0, clmBytes);  p1-=clmBytes;
		if(p0>=p1) break;
	}
	memcpy(p0, p0+clmBytes, p11-p0);
	memcpy(p11, ptmp, clmBytes);
	//Deletep(ptmp);
	return p;
}

//******************************  CContinuousCnt *******************************/
CContinuousCnt::CContinuousCnt(int cntLargest)
{
	cnt0=cntLargest0=cntLargest; isFirst=1; lost=0;
}

int CContinuousCnt::Update(int cnt)		// cnt = 0 ~ cntLargest0
{
	if(isFirst) { cnt0=cnt; isFirst=0; return 1; }
	int dcnt=cnt-cnt0; cnt0=cnt;
	if(dcnt<0) dcnt+=cntLargest0+1;
	lost += dcnt; lost--;
	return dcnt;
}

//***************************  function AlignCoarse  *********************************/
CVect3 Alignsb(const CVect3 &wmm, const CVect3 &vmm, double latitude)
{
//	if(latitude>PI/2)  latitude = asin(dot(wmm,vmm)/norm(wmm)/norm(vmm))/DEG;
	double T11, T12, T13, T21, T22, T23, T31, T32, T33;
	double cl = cos(latitude), tl = tan(latitude), nn;
	CVect3 wbib = wmm / norm(wmm),  fb = vmm / norm(vmm);
	T31 = fb.i,				 T32 = fb.j,			 	T33 = fb.k;
	T21 = wbib.i/cl-T31*tl,	 T22 = wbib.j/cl-T32*tl,	T23 = wbib.k/cl-T33*tl;		nn = sqrt(T21*T21+T22*T22+T23*T23);  T21 /= nn, T22 /= nn, T23 /= nn;
	T11 = T22*T33-T23*T32,	 T12 = T23*T31-T21*T33,		T13 = T21*T32-T22*T31;		nn = sqrt(T11*T11+T12*T12+T13*T13);  T11 /= nn, T12 /= nn, T13 /= nn;
	CMat3 Cnb(T11, T12, T13, T21, T22, T23, T31, T32, T33);
	return m2att(Cnb);
}

CVect3 Alignsb(const CVect3 &wmm, const CVect3 &vmm, const CVect3 &pos)
{
	return Alignsb(wmm, vmm, pos.i);
}

//***************************  CAligni0  *********************************/
CAligni0::CAligni0(const CVect3 &pos00, const CVect3 &vel00, int velAid0)
{
	Init(pos00, vel00, velAid0);
}

void CAligni0::Init(const CVect3 &pos00, const CVect3 &vel00, int velAid0)
{
	eth.Update(pos00);
	this->pos0 = pos00, this->vel0 = vel00, velAid = velAid0;
	tk = 0;
	t0 = t1 = 10, t2 = 0; 
	wmm = vmm = vib0 = vi0 = Pib01 = Pib02 = Pi01 = Pi02 = O31;
	qib0b = CQuat(1.0);
}

CQuat CAligni0::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, const CVect3 &vel)
{
	double nts = nSamples*ts;
	imu.Update(pwm, pvm, nSamples, ts);
	wmm = wmm + imu.phim;  vmm = vmm + imu.dvbm;
	// vtmp = qib0b * (vm + 1/2 * wm X vm)
	CVect3 vtmp = qib0b*imu.dvbm;
	// vtmp1 = qni0' * [dvn+(wnin+wnie)Xvn-gn] * ts;
	tk += nts;
	CMat3 Ci0n = pos2Cen(CVect3(eth.pos.i,eth.wie*tk,0.0));
	CVect3 vtmp1 = Ci0n*(-eth.gn*nts);
	if(velAid>0)
	{
		CVect3 dv=vel-vel0;  vel0 = vel;
		if(velAid==1)		vtmp1 += Ci0n*dv;				// for GPS vn-aided
		else if(velAid==2)	vtmp -= qib0b*dv+imu.phim*vel;	// for OD vb-aided
	}
	// Pib02 = Pib02 + vib0*ts, Pi02 = Pi02 + vi0*ts
	vib0 = vib0 + vtmp,		 vi0 = vi0 + vtmp1;
	Pib02 = Pib02 + vib0*nts, Pi02 = Pi02 + vi0*nts;
	//
	if(++t2>3*t0)
	{
		t0 = t1, Pib01 = tmpPib0, Pi01 = tmpPi0;
	}
	else if(t2>2*t0 && t1==t0)
	{
		t1 = t2, tmpPib0 = Pib02, tmpPi0 = Pi02;
	}
	//
	qib0b = qib0b*rv2q(imu.phim);
	// qnb=qni0*qiib0*qib0b
	qnbsb = a2qua(Alignsb(wmm, vmm, eth.pos.i));
	if(t2<100)
	{
		qnb0 = qnb = CQuat(1.0);
	}
	else if(t2<1000)
	{
		qnb0 = qnb = qnbsb;
	}
	else
	{
		CQuat qi0ib0 = a2qua(dv2att(Pi01, Pi02, Pib01, Pib02));
		qnb0 = (~m2qua(pos2Cen(CVect3(eth.pos.i,0.0,0.0))))*qi0ib0;
		qnb = (~m2qua(Ci0n))*qi0ib0*qib0b;
	}
	return qnb;
}



unsigned short swap16(unsigned short ui16)
{
	unsigned char *p = (unsigned char*)&ui16, c;
	c = p[0]; p[0] = p[1]; p[1] = c;
	return ui16;
}

unsigned char* swap24(unsigned char* puc3, unsigned char* pres)
{
	static unsigned char resc[3];
	if (pres == NULL) pres = resc;
	pres[0] = puc3[2]; pres[1] = puc3[1];  pres[2] = puc3[0];
	return pres;
}

unsigned int swap32(unsigned int ui32)
{
	unsigned char *p = (unsigned char*)&ui32, c;
	c = p[0]; p[0] = p[3]; p[3] = c;
	c = p[1]; p[1] = p[2]; p[2] = c;
	return ui32;
}

// unsigned long swap64(unsigned long ui64)
// {
// 	unsigned char *p = (unsigned char*)&ui64, c;
// 	c = p[0]; p[0] = p[7]; p[7] = c;
// 	c = p[1]; p[1] = p[6]; p[6] = c;
// 	c = p[2]; p[2] = p[5]; p[5] = c;
// 	c = p[3]; p[3] = p[4]; p[4] = c;
// 	return ui64;
// }

unsigned char* int24(unsigned char *pchar3, int int32)
{
	unsigned char *p = (unsigned char*)&int32;
	*pchar3++ = *p++, *pchar3++ = *p++, *pchar3 = *p;
	return pchar3-2;
}

int diffint24(const unsigned char *pc1, const unsigned char *pc0)
{
	int i1, i0;
	unsigned char *p1 = (unsigned char*)&i1, *p0 = (unsigned char*)&i0;
	*p1++ = 0, *p1++ = *pc1++, *p1++ = *pc1++, *p1 = *pc1;
	*p0++ = 0, *p0++ = *pc0++, *p0++ = *pc0++, *p0 = *pc0;
	return (i1 - i0) / 256;
}

#ifdef PSINS_psinsassert

#pragma message("  psinsassert();")

BOOL psinsassert(BOOL b)
{
	int res;

	if(b)	{
		res = 1;
	}
	else	{
		res = 0;
	}
	return res;
}

#endif

double r2dm(double r)	// rad to deg/min, eg. 1234.56 = 12deg+34.56min
{
	int sgn=1;
	if(r<0.0) { r=-r; sgn=0; }
	double deg = r/DEG;
	int ideg = (int)deg;
	double dm = ideg*100 + (deg-ideg)*60.0;
	return sgn ? dm : -dm;
}

double dm2r(double dm)
{
	int sgn=1;
	if(dm<0.0) { dm=-dm; sgn=0; }
	int ideg = (int)(dm/100);
	double r = ideg*DEG + (dm-ideg*100)*(DEG/60);
	return sgn ? r : -r;
}

BOOL logtrigger(int n, double f0)
{
	static int trigger_count=0;
	if(++trigger_count==n || f0>EPS || f0<-EPS)
	{
		trigger_count = 0;
		return TRUE;
	}
	return FALSE;
}

inline BOOL IsZero(double f, double eps)
{
	return (f<eps && f>-eps);
}

// determine the sign of 'val' with the sensitivity of 'eps'
int sign(double val, double eps)
{
	int s;

	if(val<-eps)
	{
		s = -1;
	}
	else if(val>eps)
	{
		s = 1;
	}
	else
	{
		s = 0; 
	}
	return s;
}

// set double value 'val' between range 'minVal' and 'maxVal'
double range(double val, double minVal, double maxVal)
{
	double res;

	if(val<minVal)
	{ 
		res = minVal; 
	}
	else if(val>maxVal)	
	{ 
		res = maxVal; 
	}
	else
	{ 
		res = val;
	}
	return res;
}

double  maxn(const double *pd, int n)
{
	double m=-INF; const double *pn=&pd[n];
	for(; pd<pn; pd++) {
		if(*pd>m) m=*pd;
	}
	return m;
}

double  minn(const double *pd, int n)
{
	double m=INF; const double *pn=&pd[n];
	for(; pd<pn; pd++) {
		if(*pd<m) m=*pd;
	}
	return m;
}

double norm1(const double *pd, int n)
{
	double n1=0.0; const double *pn=&pd[n];
	for(; pd<pn; pd++) {
		if(*pd>0.0) n1+=*pd; else n1-=*pd;
	}
	return n1;
}

double norm(const double *pd, int n)
{
	double n2=0.0; const double *pn=&pd[n];
	for(; pd<pn; pd++) {
		n2 += *pd*(*pd);
	}
	return sqrt(n2);
}

double normInf(const double *pd, int n)
{
	double ninf=0.0; const double *pn=&pd[n];
	for(; pd<pn; pd++) {
		if(*pd>ninf) ninf=*pd; else if(-*pd>ninf) ninf=-*pd;
	}
	return ninf;
}

double attract(double f, double th, double center)
{
	f -= center;
	if(f>-th && f<th) {
		th = f/th;  th *= th;
		f *= th;
	}
	return (f+center);
}

double polyval(const double *p, int order, double x)
{
	double y=p[0];
	for(int k=1; k<=order; k++)  y = y*x + p[k];
	return y;
}

double atan2Ex(double y, double x)
{
	double res;

	if((sign(y)==0) && (sign(x)==0))
	{
		res = 0.0;
	}
	else
	{
		res = atan2(y, x);
	}
	return res;
}

double diffYaw(double yaw, double yaw0)
{
	double dyaw = yaw-yaw0;
	if(dyaw>=PI) dyaw-=_2PI;
	else if(dyaw<=-PI) dyaw+=_2PI;
	return dyaw;
}

double MKQt(double sR, double tau)
{
	return sR*sR*2.0/tau;
}

CVect3 MKQt(const CVect3 &sR, const CVect3 &tau)
{
	return CVect3(sR.i*sR.i*2.0/tau.i, sR.j*sR.j*2.0/tau.j, sR.k*sR.k*2.0/tau.k);
}

double unixt2gpst(double ut, int leap)
{
	return fmod(ut+leap-315964800.0, 86400.0);
}

BOOL chkhdr(const char *str, const char *hdr)
{
	while(1) {
		if(*hdr=='\0') return 1;
		if(*hdr++!=*str++) return 0;
	}
}

double randn(double mu, double sigma)
{
#define NSUM 25
static double ssgm = sqrt(NSUM/12.0);
	double x = 0;
	for (int i=0; i<NSUM; i++)
	{
		x += (double)rand();
	}
	x /= RAND_MAX;
	x -= NSUM/2.0;
	x /= ssgm;		x *= sigma;
	return x+mu;
}

CVect3 randn(const CVect3 &mu, const CVect3 &sigma)
{
	return CVect3(randn(mu.i,sigma.i), randn(mu.j,sigma.j), randn(mu.k,sigma.k));
}

CVect randn(const CVect &mu, const CVect &sigma)
{
	CVect vtmp(mu.row,mu.clm);
	const double *pmu=&mu.dd[0], *pEnd=&mu.dd[mu.rc], *psigma=&sigma.dd[0];
	for(double *p=&vtmp.dd[0]; pmu<pEnd; p++,pmu++,psigma++) *p=randn(*pmu,*psigma);
	return vtmp;
}

CMat3 randn(const CMat3 &mu, const double &sigma)
{
	CMat3 mtmp;
	const double *pmu=&mu.e00, *pEnd=&mu.e22;
	for(double *p=&mtmp.e00; pmu<=pEnd; p++,pmu++) *p=randn(*pmu,sigma);
	return mtmp;
}

CMat randn(const CMat &mu, const double &sigma)
{
	CMat mtmp(mu.row,mu.clm);
	const double *pmu=&mu.dd[0], *pEnd=&mu.dd[mu.rc];
	for(double *p=&mtmp.dd[0]; pmu<pEnd; p++,pmu++) *p=randn(*pmu,sigma);
	return mtmp;
}

int* deci(int i, int *pi)
{
#define di_len 6
	static int di[di_len];   // 'i=12345' => 'di[0]=1,di[1]=2,di[2]=3,di[3]=4,di[4]=5'
	int k;
	for(k=0; k<di_len; k++) {  // decode
		di[k] = i;  i /= 10;  di[k] -= i*10;
		if(i==0) break;
	}
	for(int k1=0; k1<=k/2; k1++) {  // reverse
		int tmp=di[k1]; di[k1]=di[k-k1], di[k-k1]=tmp;
	}
	if(pi) {
		for(int j=0; j<=k; j++) pi[j]=di[j];
	}
	return di;
}

void add(CVect3 &res, const CVect3 &v1, const CVect3 &v2)
{
	res.i=v1.i+v2.i, res.j=v1.j+v2.j, res.k=v1.k+v2.k;
}

void add(CMat3 &res, const CMat3 &m1, const CMat3 &m2)
{
	res.e00=m1.e00+m2.e00, res.e01=m1.e01+m2.e01, res.e02=m1.e02+m2.e02; 
	res.e10=m1.e10+m2.e10, res.e11=m1.e11+m2.e11, res.e12=m1.e12+m2.e12; 
	res.e20=m1.e20+m2.e20, res.e21=m1.e21+m2.e21, res.e22=m1.e22+m2.e22; 
}

void add(CVect &res, const CVect &v1, const CVect &v2)
{
	const double *p1=v1.dd, *p2=v2.dd, *pEnd=&v1.dd[v1.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++,p2++) *p=*p1+*p2;
}

void add(CMat &res, const CMat &m1, const CMat &m2)
{
	const double *p1=m1.dd, *p2=m2.dd, *pEnd=&m1.dd[m1.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++,p2++) *p=*p1+*p2;
}

void subtrac(CVect3 &res, const CVect3 &v1, const CVect3 &v2)
{
	res.i=v1.i-v2.i, res.j=v1.j-v2.j, res.k=v1.k-v2.k;
}

void subtrac(CMat3 &res, const CMat3 &m1, const CMat3 &m2)
{
	res.e00=m1.e00-m2.e00, res.e01=m1.e01-m2.e01, res.e02=m1.e02-m2.e02; 
	res.e10=m1.e10-m2.e10, res.e11=m1.e11-m2.e11, res.e12=m1.e12-m2.e12; 
	res.e20=m1.e20-m2.e20, res.e21=m1.e21-m2.e21, res.e22=m1.e22-m2.e22; 
}

void subtrac(CVect &res, const CVect &v1, const CVect &v2)
{
	const double *p1=v1.dd, *p2=v2.dd, *pEnd=&v1.dd[v1.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++,p2++) *p=*p1-*p2;
}

void subtrac(CMat &res, const CMat &m1, const CMat &m2)
{
	const double *p1=m1.dd, *p2=m2.dd, *pEnd=&m1.dd[m1.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++,p2++) *p=*p1-*p2;
}

void cros(CVect3 &res, const CVect3 &v1, const CVect3 &v2)
{
	double
		i=v1.j*v2.k-v1.k*v2.j,
		j=v1.k*v2.i-v1.i*v2.k;
	res.k=v1.i*v2.j-v1.j*v2.i;
	res.i=i, res.j=j;
}

void dotmul(CVect3 &res, const CVect3 &v1, const CVect3 &v2)
{
	res.i=v1.i*v2.i, res.j=v1.j*v2.j, res.k=v1.k*v2.k;
}

void dotdiv(CVect3 &res, const CVect3 &v1, const CVect3 &v2)
{
	res.i=v1.i/v2.i, res.j=v1.j/v2.j, res.k=v1.k/v2.k;
}

void mul(CVect3 &res, const CVect3 &v, const double &f)
{
	res.i=v.i*f, res.j=v.j*f, res.k=v.k*f;
}

void mul(CMat3 &res, const CMat3 &m, const double &f)
{
	res.e00=m.e00*f, res.e01=m.e01*f, res.e02=m.e02*f;
	res.e10=m.e10*f, res.e11=m.e11*f, res.e12=m.e12*f;
	res.e20=m.e20*f, res.e21=m.e21*f, res.e22=m.e22*f;
}

void mul(CVect &res, const CVect &v, const double &f)
{
	res.row=v.row, res.clm=v.clm, res.rc=v.rc;
	const double *p1=v.dd, *pEnd=&v.dd[v.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++) *p=*p1*f;
}

void mul(CVect &res, const CMat &m, const CVect &v)
{
	res.row=res.rc=m.row, res.clm=1;
	double *p=res.dd, *pEnd=&res.dd[res.row]; const double *p1ij=m.dd, *p2End=&v.dd[v.rc];
	for(; p<pEnd; p++)
	{
		double f=0.0; const double *p2j=v.dd;
		for(; p2j<p2End; p1ij++,p2j++)	f += (*p1ij) * (*p2j);
		*p = f;
	}
}

void mul(CMat &res, const CMat &m, const double &f)
{
	res.row=m.row, res.clm=m.clm, res.rc=m.rc;
	const double *p1=m.dd, *pEnd=&m.dd[m.rc];
	for(double *p=res.dd; p1<pEnd; p++,p1++) *p=*p1*f;
}

void mul(CVect3 &res, const CMat3 &m, const CVect3 &v)
{
	res.i=m.e00*v.i+m.e01*v.j+m.e02*v.k, res.j=m.e10*v.i+m.e11*v.j+m.e12*v.k, res.k=m.e20*v.i+m.e21*v.j+m.e22*v.k;
}

void mul(CVect3 &res, const CVect3 &v, const CMat3 &m)
{
	res.i=m.e00*v.i+m.e10*v.j+m.e20*v.k, res.j=m.e01*v.i+m.e11*v.j+m.e21*v.k, res.k=m.e02*v.i+m.e12*v.j+m.e22*v.k;
}

void mul(CMat3 &res, const CMat3 &m1, const CMat3 &m2)
{
	if(&res==&m1||&res==&m2) {
		res=m1*m2;
	}
	else {
		res.e00 = m1.e00*m2.e00 + m1.e01*m2.e10 + m1.e02*m2.e20;
		res.e01 = m1.e00*m2.e01 + m1.e01*m2.e11 + m1.e02*m2.e21;
		res.e02 = m1.e00*m2.e02 + m1.e01*m2.e12 + m1.e02*m2.e22;
		res.e10 = m1.e10*m2.e00 + m1.e11*m2.e10 + m1.e12*m2.e20;
		res.e11 = m1.e10*m2.e01 + m1.e11*m2.e11 + m1.e12*m2.e21;
		res.e12 = m1.e10*m2.e02 + m1.e11*m2.e12 + m1.e12*m2.e22;
		res.e20 = m1.e20*m2.e00 + m1.e21*m2.e10 + m1.e22*m2.e20;
		res.e21 = m1.e20*m2.e01 + m1.e21*m2.e11 + m1.e22*m2.e21;
		res.e22 = m1.e20*m2.e02 + m1.e21*m2.e12 + m1.e22*m2.e22;
	}
}

void mul(CQuat &res, const CQuat &q1, const CQuat &q2)
{
	double 
		qq0 = q1.q0*q2.q0 - q1.q1*q2.q1 - q1.q2*q2.q2 - q1.q3*q2.q3,
		qq1 = q1.q0*q2.q1 + q1.q1*q2.q0 + q1.q2*q2.q3 - q1.q3*q2.q2,
		qq2 = q1.q0*q2.q2 + q1.q2*q2.q0 + q1.q3*q2.q1 - q1.q1*q2.q3;
	 res.q3 = q1.q0*q2.q3 + q1.q3*q2.q0 + q1.q1*q2.q2 - q1.q2*q2.q1;
	res.q0=qq0, res.q1=qq1, res.q2=qq2;
}

void mul(CVect3 &res, const CQuat &q, const CVect3 &v)
{
	double
		qq0 =         - q.q1*v.i - q.q2*v.j - q.q3*v.k,
		qq1 = q.q0*v.i           + q.q2*v.k - q.q3*v.j,
		qq2 = q.q0*v.j           + q.q3*v.i - q.q1*v.k,
		qq3 = q.q0*v.k           + q.q1*v.j - q.q2*v.i;
	res.i = -qq0*q.q1 + qq1*q.q0 - qq2*q.q3 + qq3*q.q2;
	res.j = -qq0*q.q2 + qq2*q.q0 - qq3*q.q1 + qq1*q.q3;
	res.k = -qq0*q.q3 + qq3*q.q0 - qq1*q.q2 + qq2*q.q1;
}

void _TT(CMat3 &mT, const CMat3 &m)
{
	mT.e00=m.e00, mT.e10=m.e01, mT.e20=m.e02;
	mT.e01=m.e10, mT.e11=m.e11, mT.e21=m.e12;
	mT.e02=m.e20, mT.e12=m.e21, mT.e22=m.e22;
}

void rv2q(CQuat &q, const CVect3 &rv)
{
	double n2 = rv.i*rv.i+rv.j*rv.j+rv.k*rv.k, f;
	if(n2<(PI/180.0*PI/180.0))	// 0.017^2 
	{
		double n4=n2*n2;
		q.q0 = 1.0 - n2*(1.0/rvF2) + n4*(1.0/rvF4);
		f = 0.5 - n2*(1.0/rvF3) + n4*(1.0/rvF5);
	}
	else
	{
		double n_2 = sqrt(n2)/2.0;
		q.q0 = cos(n_2);
		f = sin(n_2)/n_2*0.5;
	}
	q.q1=f*rv.i, q.q2=f*rv.j, q.q3=f*rv.k;
}

void qdelphi(CQuat &q, const CVect3 &phi)
{
	CQuat qtmp;
	rv2q(qtmp, phi);
	double
		  q0 = qtmp.q0*q.q0 - qtmp.q1*q.q1 - qtmp.q2*q.q2 - qtmp.q3*q.q3,
		  q1 = qtmp.q0*q.q1 + qtmp.q1*q.q0 + qtmp.q2*q.q3 - qtmp.q3*q.q2,
		  q2 = qtmp.q0*q.q2 + qtmp.q2*q.q0 + qtmp.q3*q.q1 - qtmp.q1*q.q3;
		q.q3 = qtmp.q0*q.q3 + qtmp.q3*q.q0 + qtmp.q1*q.q2 - qtmp.q2*q.q1;
	q.q0=q0; q.q1=q1; q.q2=q2;
}

void q2mat(CMat3 &Cnb, const CQuat &qnb)
{
	double	q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
			q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
			q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
			q44 = qnb.q3*qnb.q3;
    Cnb.e00 = q11+q22-q33-q44,  Cnb.e01 = 2*(q23-q14),     Cnb.e02 = 2*(q24+q13),
	Cnb.e10 = 2*(q23+q14),      Cnb.e11 = q11-q22+q33-q44, Cnb.e12 = 2*(q34-q12),
	Cnb.e20 = 2*(q24-q13),      Cnb.e21 = 2*(q34+q12),     Cnb.e22 = q11-q22-q33+q44;
}

void m2att(CVect3 &att, const CMat3 &Cnb)
{
//	att.i = asinEx(Cnb.e21);
	float e21=(float)Cnb.e21;
	if(e21<-1.0f) e21=1.0; else if(e21>1.0f) e21=1.0;  att.i = asinf(e21);
//	att.j = atan2Ex(-Cnb.e20, Cnb.e22);
	float e20=(float)Cnb.e20, e22=(float)Cnb.e22;
	att.j = ((double)e20>1.0e-10||(double)e20<-1.0e-10 || (double)e22>1.0e-10||(double)e22<-1.0e-10) ? (double)atan2f(-e20, e22) : 0.0;
//	att.k = atan2Ex(-Cnb.e01, Cnb.e11);
	float e01=(float)Cnb.e01, e11=(float)Cnb.e11;
	att.k = ((double)e01>1.0e-10||(double)e01<-1.0e-10 || (double)e11>1.0e-10||(double)e11<-1.0e-10) ? (double)atan2f(-e01, e11) : 0.0;
}

void AXbt(CVect3 &res, const CMat3 &A, const CVect3 &X, const CVect3 &b, const double &t)
{
	double  i = A.e00*X.i+A.e01*X.j+A.e02*X.k + b.i*t,
		    j = A.e10*X.i+A.e11*X.j+A.e12*X.k + b.j*t;
		res.k = A.e20*X.i+A.e21*X.j+A.e22*X.k + b.k*t;
	res.i=i, res.j=j;
}

void sizedisp(int i)
{
#define ClsSize(xxx)	printf("\t"#xxx":%10ld\n", sizeof(xxx))
	printf("Size of each class (in bytes):\n");
	ClsSize(CVect);
	ClsSize(CMat);
	ClsSize(CMaxMinn);
	ClsSize(CAligni0);
	ClsSize(CSINS);
	ClsSize(CKalman);
	ClsSize(CSINSTDKF);
	ClsSize(CSINSGNSS);
	ClsSize(CRAvar);
	if(i>0) printf("\tclassXXX: %10d\n", i);
	exit(0);
}
