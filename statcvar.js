let igammaepsilon = 0.000000000000001;
let igammabignumber = 4503599627370496.0;     
let igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
let maxrealnumber=Math.pow(10,300);
let minrealnumber=Math.pow(10,-300);
let s2pi=Math.sqrt(2*Math.PI);

//Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

function invnormaldistribution(y0)  {

let expm2 = 0.13533528323661269189;
let s2pi = 2.50662827463100050242;
            let result = 0;
            let x = 0;
            let y = 0;
            let z = 0;
            let y2 = 0;
            let x0 = 0;
            let x1 = 0;
            let code = 0;
            let p0 = 0;
            let q0 = 0;
            let p1 = 0;
            let q1 = 0;
            let p2 = 0;
            let q2 = 0;
            let zz;
            
            if(y0<=0) {
                result = -maxrealnumber;
                return result;
            }
            if(y0>=1) {
                result = maxrealnumber;
                return result;
            }
            code = 1;
            y = y0;
            if(y>1.0-expm2) {
                y = 1.0-y;
                code = 0;
                }
            if(y>expm2) {
                y = y-0.5;
                y2 = y*y;
                p0 = -59.9633501014107895267;
                p0 = 98.0010754185999661536+y2*p0;
                p0 = -56.6762857469070293439+y2*p0;
                p0 = 13.9312609387279679503+y2*p0;
                p0 = -1.23916583867381258016+y2*p0;
                q0 = 1;
                q0 = 1.95448858338141759834+y2*q0;
                q0 = 4.67627912898881538453+y2*q0;
                q0 = 86.3602421390890590575+y2*q0;
                q0 = -225.462687854119370527+y2*q0;
                q0 = 200.260212380060660359+y2*q0;
                q0 = -82.0372256168333339912+y2*q0;
                q0 = 15.9056225126211695515+y2*q0;
                q0 = -1.18331621121330003142+y2*q0;
                x = y+y*y2*p0/q0;
                x = x*s2pi;
                result = x;
                return result;
            }
            zz=Math.log(y);
            x = Math.sqrt(-(2.0*zz));
            x0 = x-Math.log(x)/x;
            z = 1.0/x;
            if(x<8.0) {
                p1 = 4.05544892305962419923;
                p1 = 31.5251094599893866154+z*p1;
                p1 = 57.1628192246421288162+z*p1;
                p1 = 44.0805073893200834700+z*p1;
                p1 = 14.6849561928858024014+z*p1;
                p1 = 2.18663306850790267539+z*p1;
                p1 = -(1.40256079171354495875*0.1)+z*p1;
                p1 = -(3.50424626827848203418*0.01)+z*p1;
                p1 = -(8.57456785154685413611*0.0001)+z*p1;
                q1 = 1;
                q1 = 15.7799883256466749731+z*q1;
                q1 = 45.3907635128879210584+z*q1;
                q1 = 41.3172038254672030440+z*q1;
                q1 = 15.0425385692907503408+z*q1;
                q1 = 2.50464946208309415979+z*q1;
                q1 = -(1.42182922854787788574*0.1)+z*q1;
                q1 = -(3.80806407691578277194*0.01)+z*q1;
                q1 = -(9.33259480895457427372*0.0001)+z*q1;
                x1 = z*p1/q1;
            }
            else  {
                p2 = 3.23774891776946035970;
                p2 = 6.91522889068984211695+z*p2;
                p2 = 3.93881025292474443415+z*p2;
                p2 = 1.33303460815807542389+z*p2;
                p2 = 2.01485389549179081538*0.1+z*p2;
                p2 = 1.23716634817820021358*0.01+z*p2;
                p2 = 3.01581553508235416007*0.0001+z*p2;
                p2 = 2.65806974686737550832*0.000001+z*p2;
                p2 = 6.23974539184983293730*0.000000001+z*p2;
                q2 = 1;
                q2 = 6.02427039364742014255+z*q2;
                q2 = 3.67983563856160859403+z*q2;
                q2 = 1.37702099489081330271+z*q2;
                q2 = 2.16236993594496635890*0.1+z*q2;
                q2 = 1.34204006088543189037*0.01+z*q2;
                q2 = 3.28014464682127739104*0.0001+z*q2;
                q2 = 2.89247864745380683936*0.000001+z*q2;
                q2 = 6.79019408009981274425*0.000000001+z*q2;
                x1 = z*p2/q2;
            }
            x = x0-x1;
            if( code!=0 ) x = -x;
            result = x;
            return result;
          }
//************************************************************************
 function errorfunctionc(x) {
            let result = 0;
            let p = 0;
            let q = 0;
   
            if(x<0) {
                result = 2-errorfunctionc(-x);
                return result;
            }
            if(x<0.5) {
                result = 1.0-errorfunction(x);
                return result;
            }
            if(x>=10) {
                result = 0;
                return result;
            }
            p = 0.0;
            p = 0.5641877825507397413087057563+x*p;
            p = 9.675807882987265400604202961+x*p;
            p = 77.08161730368428609781633646+x*p;
            p = 368.5196154710010637133875746+x*p;
            p = 1143.262070703886173606073338+x*p;
            p = 2320.439590251635247384768711+x*p;
            p = 2898.0293292167655611275846+x*p;
            p = 1826.3348842295112592168999+x*p;
            q = 1.0;
            q = 17.14980943627607849376131193+x*q;
            q = 137.1255960500622202878443578+x*q;
            q = 661.7361207107653469211984771+x*q;
            q = 2094.384367789539593790281779+x*q;
            q = 4429.612803883682726711528526+x*q;
            q = 6089.5424232724435504633068+x*q;
            q = 4958.82756472114071495438422+x*q;
            q = 1826.3348842295112595576438+x*q;
            result = Math.exp(-x*x)*p/q;
            return result;
        }

//*****************************************************************
 function errorfunction(x) {
            let result = 0;
            let xsq = 0;
            let s = 0;
            let p = 0;
            let q = 0;

            s = Math.sign(x);
            x = Math.abs(x);
            if(x<0.5) {
                xsq = x*x;
                p = 0.007547728033418631287834;
                p = -0.288805137207594084924010+xsq*p;
                p = 14.3383842191748205576712+xsq*p;
                p = 38.0140318123903008244444+xsq*p;
                p = 3017.82788536507577809226+xsq*p;
                p = 7404.07142710151470082064+xsq*p;
                p = 80437.3630960840172832162+xsq*p;
                q = 0.0;
                q = 1.0+xsq*q;
                q = 38.0190713951939403753468+xsq*q;
                q = 658.070155459240506326937+xsq*q;
                q = 6379.60017324428279487120+xsq*q;
                q = 34216.5257924628539769006+xsq*q;
                q = 80437.3630960840172826266+xsq*q;
                result = s*1.1283791670955125738961589031*x*p/q;
                return result;
            }
            if(x>=10) {
                result = s;
                return result;
            }
            result = s*(1-errorfunctionc(x));
            return result;
        }

//**************************************************************************
 function normaldistribution(x) {
            let result=0;
            result=0.5*(errorfunction(x/1.41421356237309504880)+1);
            return result;
        }
//*******************************************************************
function gammastirf(x) {
 
let stir,w,y,v,result; 
w = 1./x;
stir = 7.87311395793093628397*0.0001;
stir = -2.29549961613378126380*0.0001+w*stir;
stir = -2.68132617805781232825*0.001+w*stir;
stir = 3.47222221605458667310*0.001+w*stir;
stir = 0.0833333333333482257126+w*stir;
w = 1.+w*stir;y = Math.exp(x);

if (x>143.01608)  {
   v = Math.pow(x, 0.5*x-0.25);
   y = v*v/y;
   }
 else  {
   y = Math.pow(x, x-0.5)/y;
   }
   result =2.50662827463100050242*y*w;
   return result;
}
//*****************************************************************************
function gamma(x) {
  let sgngam=1;
  let q=Math.abs(x);
  let i,p,z,result,pp,qq;

if(q>33) {
  if(x<0) {
    p=Math.floor(q);
    i=Math.round(p);
    if(i%2==0) sgngam=-1;
    z=q-p;
    if(z>0.5) {
      p++;
      z=q-p;
    }
    z=q*Math.sin(Math.PI*z);
    z=Math.abs(z);
    z=Math.PI/(z*gammastirf(q));
  } else {
    z=gammastirf(x);
  }
  result=sgngam*z;
  return result;
}

z=1;
while(x>=3) {
  x--;
  z*=x;
}
while(x<0) {
  if(x>-0.000000001) {
    result=z/((1+0.577215664901533*x)*x);
    return result;
  }
  z/=x;
  x++;
}
while(x<2) {
  if(x<0.000000001) {
    result = z/((1 + 0.577215664901533*x)*x);
    return result;
  }
  z/=x;
  x++;
}

if(x==2) {
  result=z;
  return result;
}
    x = x - 2;
    pp = 1.60119522476752*0.0001;
    pp = 1.19135147006586*0.001 + x * pp;
    pp = 1.04213797561762*0.01 + x * pp;
    pp = 4.76367800457137*0.01 + x * pp;
    pp = 0.207448227648436 + x * pp;
    pp = 0.494214826801497 + x * pp;
    pp = 1 + x*pp;
    qq = -2.3158187332412*0.00001;
    qq = 5.39605580493303*0.0001 + x*qq;
    qq = -4.45641913851797*0.001 + x*qq;
    qq = 0.011813978522206 + x*qq;
    qq = 3.58236398605499*0.01 + x*qq;
    qq = -0.234591795718243 + x*qq;
    qq = 7.14304917030273*0.01 + x*qq;
    qq = 1 + x*qq;
    result = z * pp / qq;
    return result;
}

//*****************************************************************************
function invchisquaredistribution(v,y) {
 let result;
 result = 2*invincompletegammac(0.5*v,y);
 return(result);
}
//************************************************************************
function invincompletegammac(a,y0) {
            
            let i,result,dir;
            let tmp=0;     
            let x0 = igammabignumber;
            let yl = 0;
            let x1 = 0;
            let yh = 1;
            let dithresh = 5*igammaepsilon;
            let d = 1/(9*a);
            let y = 1-d-invnormaldistribution(y0)*Math.sqrt(d);
            let x = a*y*y*y;
            let lgm = lngamma(a,tmp);
            i = 0;
            while( i<10 )
            {
                if(x>x0 || x<x)
                {
                    d = 0.0625;
                    break;
                }
                y = incompletegammac(a, x);
                if(y<yl || y>yh)
                {
                    d = 0.0625;
                    break;
                }
                if(y<y0)
                {
                    x0 = x;
                    yl = y;
                }
                else
                {
                    x1 = x;
                    yh = y;
                }
                d = (a-1)*Math.log(x)-x-lgm;
                if (d<-709.78271289338399)
                {
                    d = 0.0625;
                    break;
                }
                d = -Math.exp(d);
                d = (y-y0)/d;
                if(Math.abs(d/x)<igammaepsilon)
                {
                    result = x;
                    return result;
                }
                x = x-d;
                i = i+1;
            }
            if (x0==igammabignumber) {
                if (x<=0) x = 1;
                 while(x0==igammabignumber) {
                    x = (1+d)*x;
                    y = incompletegammac(a, x);
                    if (y<y0)   {
                        x0 = x;
                        yl = y;
                        break;
                    }
                    d = d+d;
                }
            }
            d = 0.5;
            dir = 0;
            i = 0;
            while (i<400) {
                x = x1+d*(x0-x1);
                y = incompletegammac(a, x);
                lgm = (x0-x1)/(x1+x0);
                if (Math.abs(lgm)<dithresh) break;
                lgm = (y-y0)/y0;
                if (Math.abs(lgm)<dithresh) break;
                if (x<=0)  break;
                if(y>=y0) {
                    x1 = x;
                    yh = y;
                    if (dir<0) {
                        dir = 0;d = 0.5;
                    }
                    else
                    {
                        if (dir>1) {
                            d = 0.5*d+0.5;
                        }
                        else  {
                            d = (y0-yl)/(yh-yl);
                        }
                    }
                    dir = dir+1;
                }
                else
                {
                    x0 = x;
                    yl = y;
                    if (dir>0)
                    {
                        dir = 0;
                        d = 0.5;
                    }
                    else
                    {
                        if (dir<-1)
                        {
                            d = 0.5*d;
                        }
                        else
                        {
                            d = (y0-yl)/(yh-yl);
                        }
                    }
                    dir = dir-1;
                }
                i = i+1;
            }
            result = x;
            return result;
        }

//**************************************************************************
function incompletegamma(a,x) {

  let tmp=0;
  let ax,r,c,ans,result;
 
if(x<=0 || a<=0) {
  result=0;
  return result;
}
if(x>1 && x>a) {
   result=1-incompletegammac(a, x);
   return result;
}
ax=a*Math.log(x)-x-lngamma(a,tmp);
if(ax<-709.782712893384) {
  result=0;
  return result;
}
ax=Math.exp(ax);
r=a;c=1;ans=1;
do {
  r++;c*=x/r;ans+=c;
} while((c/ans)>igammaepsilon);
result=ans*ax/a;
return result;
}
//*****************************************************************************
 function incompletegammac(a,x) {
     
        let tmp=0;
        let result,ax,y,z,c,ans,pkm2,qkmw,pkm1,qkm1,yc,pk,qk,r,t;

        if(x<=0 || a<=0) {
                result = 1;
                return result;
            }
            if( x<1 || x<a)
            {
                result = 1-incompletegamma(a, x);
                return result;
            }
            ax = a*Math.log(x)-x-lngamma(a,tmp);
            if(ax<-709.78271289338399)
            {
                result = 0;
                return result;
            }
            ax = Math.exp(ax);
            y = 1-a;
            z = x+y+1;
            c = 0;
            pkm2 = 1;
            qkm2 = x;
            pkm1 = x+1;
            qkm1 = z*x;
            ans = pkm1/qkm1;
            do
            {
                c = c+1;
                y = y+1;
                z = z+2;
                yc = y*c;
                pk = pkm1*z-pkm2*yc;
                qk = qkm1*z-qkm2*yc;
                if (qk!=0)
                {
                    r = pk/qk;
                    t = Math.abs((ans-r)/r);
                    ans = r;
                }
                else
                {
                    t = 1;
                }
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if(Math.abs(pk)>igammabignumber)
                {
                    pkm2 = pkm2*igammabignumberinv;
                    pkm1 = pkm1*igammabignumberinv;
                    qkm2 = qkm2*igammabignumberinv;
                    qkm1 = qkm1*igammabignumberinv;
                }
            } while(t>igammaepsilon);
            result = ans*ax;
            return result;
        }
//*************************************************************
  function chisquaredistribution(v,x) {
            let result;

            result = incompletegamma(v/2.0, x/2.0);
            return result;
        }
//************************************************************
function lngamma(x,sgngam) {

     sgngam = 1;
     let q,p,i,z,w,u,result,a,b,c;
    
    if(x<-34) {
      q=-x;
      w=lngamma(q,tmp);
      p=Math.floor(q);
      i=Math.round(p);
      if(i%2==0) {
        sgngam=-1;
      } else {
        sgngam=1;
      }
      z=q-p;
      if(z>0.5) {
        p++;
        z=p-q;
      }
      z=q*Math.sin(Math.PI*z);
       result=loqpi-Math.log(z)-w;
       return result;
    }
    if(x<13) {
      z=1;
      p=0;
      u=x;
      while(u>=3) {
        p--;
        u=x+p;
        z*=u;
      }
      while(u<2) {
        z/=u;
        p++;
        u=x+p;
      }
      if(z<0) {
        sgngam=-1;
        z=-z;
      } else {
        sgngam=1;
      }
      if(u==2) {
        result=Math.log(z);
        return result;
      }
        p = p - 2;
        x = x + p;
        b = -1378.25152569121;
        b = -38801.6315134638 + x * b;
        b = -331612.992738871 + x * b;
        b = -1162370.97492762 + x * b;
        b = -1721737.0082084 + x * b;
        b = -853555.664245765 + x * b;
        c = 1;
        c = -351.815701436523 + x * c;
        c = -17064.2106651881 + x * c;
        c = -220528.590553854 + x * c;
        c = -1139334.44367983 + x * c;
        c = -2532523.07177583 + x * c;
        c = -2018891.41433533 + x * c;
        p = x * b / c;
        result = Math.log(z) + p;
        return result;
    }
    q=(x-0.5)*Math.log(x)-x+ls2pi;
    if(x>100000000) {
      result=q;
      return result;
    }
    p=1/(x*x);
    if(x>=1000) {
      q+=((7.93650793650794 * 0.0001 * p - 2.77777777777778 * 0.001) * p + 0.0833333333333333) / x;
    } else {
        a = 8.11614167470508 * 0.0001;
        a = -(5.95061904284301 * 0.0001) + p * a;
        a = 7.93650340457717 * 0.0001 + p * a;
        a = -(2.777777777301 * 0.001) + p * a;
        a = 8.33333333333332 * 0.01 + p * a;
        q = q + a / x;
    }
    result=q;
    return result;
}
//******************************************************
function sf35r(x) {
      let s,t,z;
      t =Math.abs(x);
      s=0;
      if(t<6.5) {
      s=Math.exp(-t*t)*((((((0.56419*t+6.802899)*t+38.71143)*t+131.1266)*t+278.5978)*t+355.969)*t+224.1828)/(((((((t+12.05784)*t+69.11384)*t+238.4503)*t+527.5538)*t+741.5214)*t+608.9322)*t+224.1828);
      }      
      z=s; 
      if(x<0) z=2.0-s;
      return(z);
}

//*****************************************************

function sf49r(x) {
 return(0.5*sf35r(-0.7071067*x));
}
//*******************************************************

function sf53r(y,z,eps) {

 let ep1, t, b, a, ta;
 let c, hsqb, expov,asq,a4;
 let b4,a4b4, ahsqb, ab4;
 let f, sum, gt, g1, ber;
 let ter, d1, d2, aeps, d;
 let bexp, sys076, sys017;

 sys076=88.72283;sys017=0.1591549;expov=88.72283;c=0.1591549;ep1=eps;
 if(eps==0) ep1=0.000001;
 t=0;b=Math.abs(y);a=Math.abs(z);
 if(a==0) return t;
 ta=Math.atan(a);
 

 if(a*b>4) {
  t=sf49r(b);
  t=c*(ta+Math.atan(1/a))-0.5*(t-0.5);
  if(z<0) t=-t;
  return(t);
 }


   hsqb=0.5*b*b;
   if(hsqb>expov) return(t);
   bexp=Math.exp(-hsqb);
   asq=a*a;a4=asq*asq;b4=hsqb*hsqb;a4b4=a4*b4;ahsqb=a*hsqb;ab4=a*b4*0.5;f=1;sum=0;g=3;

//***************************************
 do {
    g1=g;ber=0;ter=ab4;
 
   do {
     ber=ber+ter;
     ter=ter*hsqb/g1;
     g1=g1+1;
    } while(ter>(ber*ep1));

   d1=(ber+ahsqb)/f;
   d2=ber*asq/(f+2);
   d=d1-d2;
   sum=sum+d;
   t=ta-sum*bexp;
   aeps=ep1*t;
   ahsqb=ahsqb*a4b4/((g-1)*g);
   ab4=ab4*a4b4/((g+1)*g);
   f=f+4;
   g=g+2;
  } while((d2*bexp)>=aeps);
//********************************

   t=t*c;
   if(z<0) t=-t;
   return(t);
}

//**************************************************************     
function sf54r(x,d,idf) {

/*
Parameters:
   Input: x, the argument.
   idf, the number of degrees of freedom.
   d the noncentrality parameter.
   Output: zsf54r-lower tail of the noncentral t distribution.
*/

  let l,i1,idfm2,ss1,ss2,df,tval,a,b,sb,da;
  let dsb,dasb,p1,f2,a1,a2,a3,f1,sum,az,fz,fkm1,c,p,zsf54r; 

  c=0.1591549;a2=0.1591549;a3=2.506628;a1=0.3989423;ss1=0.00001;
  ss2=0.7071068;
  if(idf<=0) return(3.4e+38);
  df=parseFloat(idf);tval=x;
  
  i1=idf-2*Math.round((0.5*idf));
  i1=Math.abs(i1);
  
  a=tval/Math.sqrt(df);b=df/(df+tval*tval);sb=Math.sqrt(b);da=d*a;
  dsb=d*sb;dasb=a*dsb;p1=sf49r(dasb);
  f2=a*sb*Math.exp(-0.5*dsb*dsb)*p1*a1;
  f1=b*(da*f2+a*a2*Math.exp(-0.5*d*d));
 
  sum=0;
  
  if(idf!=1) {
  if(i1>0) {
    sum=f1;
  }
  else {
  sum=f2;
  }
  if(idf>=4) {
    idfm2=idf-2;az=1;fz=2;
    for(l=2;l<=idfm2;l+=2) {
      fkm1=fz-1;f2=b*(da*az*f1+f2)*fkm1/fz;
      az=1/(az*fkm1);
      f1=b*(da*az*f2+f1)*fz/(fz+1);
      if(i1<=0) {
       sum=sum+f2;
      }
      else {
        sum=sum+f1;
      }
       az=1/(az*fz);
       fz=fz+2;
    }
  }

  }
     
       if(i1>0) {
         p1=0.5*sf35r(ss2*dsb);
         p=sf53r(dsb,a,ss1);
         zsf54r=p1+2*(p+sum);
        }
        else {
         p1=0.5*sf35r(ss2*d);
         zsf54r=p1+sum*a3;
       }
       if(zsf54r<0) zsf54r=0;

  return(zsf54r);      
}
//****************************************************************************
function prncst(st,idf,d) {

/*
   Purpose:
   PRNCST computes the lower tail of noncentral T distribution.
   Licensing:
   This code is distributed under the GNU LGPL license. 
   Modified:
   26 January 2008
   Author:
   Original FORTRAN77 version by BE Cooper.
   C++ version by John Burkardt.
   Reference:
   BE Cooper,
   Algorithm AS 5:
   The Integral of the Non-Central T-Distribution,
   Applied Statistics,Volume 17, Number 2, 1968, page 193.
   Parameters:
   Input: st, the argument.
   idf, the number of degrees of freedom.
   d the noncentrality parameter.
   Output: valuex -lower tail of the noncentral t distribution.
   Local, g1=1.0/ Math.sqrt(2.0 * pi)
   g2=1.0 / (2.0 * pi)
   g3= Math.sqrt(2.0 * pi)
*/

  let a;
  let ak;
  let b;
  let da;
  let drb;
  let emin = 12.5;
  let f;
  let fk;
  let fkm1;
  let fmkm1;
  let fmkm2;
  let g1 = 0.3989422804;
  let g2 = 0.1591549431;
  let g3 = 2.5066282746;
  let ioe;
  let k;
  let rb;
  let sum;
  let valuex;

  f =parseFloat(idf);
//  For very large IDF, use the normal approximation.
  if ( 100 < idf ) {
    ifault = 1;
    a = d*Math.sqrt(0.5*f)*gamma(0.5*(f-1.0))/gamma(0.5*f);
    z=(st-a)/Math.sqrt(f*(1.+d*d)/(f-2.0)-a*a);
    valuex = normaldistribution(z);
    return valuex;
  }

  ifault = 0;
  ioe =idf % 2;
  a = st / Math.sqrt (f);
  b = f / ( f + st * st );
  rb = Math.sqrt ( b );
  da = d * a;
  drb = d * rb;
  if ( idf == 1 ) {
    valuex=1-normaldistribution(drb)+2.0*tfn(drb,a );
    return valuex;
  }
  sum = 0.0;
  if ( Math.abs ( drb ) < emin )  {
   fmkm2=a*rb*Math.exp(-0.5*drb*drb )*normaldistribution(a*drb)*g1;
}
  else {
    fmkm2 = 0.0;
  }
  fmkm1 = b * da * fmkm2;
  if ( Math.abs ( d ) < emin ) {
    fmkm1 = fmkm1 + b * a * g2 * Math.exp ( - 0.5 * d * d );
  }

  if ( ioe == 0 ) {
    sum = fmkm2;
  }
  else {
    sum = fmkm1;
  }

  ak = 1.0;
  fk = 2.0;

  for ( k = 2; k <= idf - 2; k = k + 2 )  {
    fkm1 = fk - 1.0;
    fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk;
    ak = 1.0 / ( ak * fkm1 );
    fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );
    if ( ioe == 0 ) {
      sum = sum + fmkm2;
    }
    else {
      sum = sum + fmkm1;
    }
    ak = 1.0 / ( ak * fk );
    fk = fk + 2.0;
  }
  if ( ioe == 0 ) {
  valuex =1-normaldistribution(d)+sum*g3;
  }
  else {
  valuex =1-normaldistribution(drb)+2.0*(sum+tfn(drb,a));
   }
  return valuex;
}
//****************************************************************************
function tfn(x,fx) {

/*
  Purpose:
  tfn calculates the T-function of Owen.
  Licensing:
  This code is distributed under the GNU LGPL license. 
  Modified:16 January 2008
  Author:
  Original FORTRAN77 version by JC Young, Christoph Minder.
  C++ version by John Burkardt.
  Reference:
  MA Porter, DJ Winstanley,
  Remark AS R30:
  A Remark on Algorithm AS76:
  An Integral Useful in Calculating Noncentral T and Bivariate
  Normal Probabilities,
  Applied Statistics,Volume 28, Number 1, 1979, page 113.
  JC Young, Christoph Minder,
  Algorithm AS 76: 
  An Algorithm Useful in Calculating Non-Central T and 
  Bivariate Normal Distributions,
  Applied Statistics,Volume 23, Number 3, 1974, pages 455-457.
  Parameters:
  Input: x, fx the parameters of the function.
  Output: valuex of the t-function.
*/

  ng=5;

  let fxs;
  let i;
  let r=[0.1477621,0.1346334,0.1095432,0.0747257,0.0333357];
  let r1;
  let r2;
  let rt;
  let tp = 0.159155;
  let tv1 = 1.0E-35;
  let tv2 = 15.0;
  let tv3 = 15.0;
  let tv4 = 1.0E-05;
  let u=[0.0744372,0.2166977,0.3397048,0.4325317,0.4869533];
  let valuex;
  let x1;
  let x2;
  let xs;

//  Test for X near zero.

  if ( Math.abs ( x ) < tv1 ) {
    valuex = tp * Math.atan ( fx );
    return valuex;
  }

//  Test for large values of abs(x).

  if ( tv2 < Math.abs ( x ) ) {
    valuex = 0.0;
    return valuex;
  }

//  Test for fx near zero.

  if ( Math.abs ( fx ) < tv1 ) {
    valuex = 0.0;
    return valuex;
  }

//  Test whether abs(fx) is so large that it must be truncated.

  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
//  Computation of truncation point by Newton iteration.
  if ( tv3 <= Math.log ( 1.0 + fxs ) - xs * fxs ){
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;
    for ( ; ; )  {
      rt = fxs + 1.0;
      x2 = x1 + ( xs * fxs + tv3 - Math.log ( rt ) )/ ( 2.0 * x1 * ( 1.0 / rt - xs ) );
      fxs = x2 * x2;
      if ( Math.abs ( x2 - x1 ) < tv4 ) {
        break;
      }
      x1 = x2;
    }
  }
//  Gaussian quadrature.
  rt = 0.0;
  for ( i = 0; i < ng; i++ ){
    r1 = 1.0 + fxs * Math.pow ( 0.5 + u[i], 2 );
    r2 = 1.0 + fxs * Math.pow ( 0.5 - u[i], 2 );
    rt = rt + r[i] * ( Math.exp ( xs * r1 ) / r1 + Math.exp ( xs * r2 ) / r2 );
  }
  valuex = rt * x2 * tp;
  return valuex;
}
//***********Вейбулловские порядковые статистики*************************************
//Copyright (c) 1994-2021, Levon Agamirov

function orderw(n,pr,ps) {

/*
 Computes David-Johnson approximation for mean and covariance between rth
 and sth order statistics from the Weibull dist. for a sample size n.
 pr=r/(n+1);ps=s/(n+1);
 vorder[0] - mean;
 vorder[1] - covariance between rth and sth order statistics 
  if r=s then vorder[1] - variance
*/

  let qr,qs,xr,xs,dr,ds;
  let xr1,xr2,xr3,xr4,xr5,xr55,xr6;
  let xs1,xs2,xs3,xs4,xs5,xs55,xs6;
  let z1,z2,z3,z4,z5,z6,z7;
  let a1,b1,c1,d1;
  let er,crs;


 qr=1-pr; qs=1-ps;
xr=Math.log(Math.log(1/(1-pr)));
xs=Math.log(Math.log(1/(1-ps)));
xr1=1/(Math.log(1/(1-pr))*(1-pr));
xr2=xr1*(1/(1-pr)-xr1);
xr3=xr2*xr2/xr1+xr1*(1/Math.pow(1-pr,2)-xr2);
xr4=(3*xr1*xr2*xr3-2*Math.pow(xr2,3))/Math.pow(xr1,2)+xr1*(2/Math.pow(1-pr,3)-xr3);
xr55=(-12*xr1*Math.pow(xr2,2)*xr3+3*Math.pow(xr1,2)*Math.pow(xr3,2)+4*Math.pow(xr1,2)*xr2*xr4+6*Math.pow(xr2,4));
xr5=xr55/Math.pow(xr1,3)+xr1*(6/Math.pow(1-pr,4)-xr4);
a1=-12*Math.pow(xr2,3)*xr3-12*xr1*(2*xr2*Math.pow(xr3,2)+Math.pow(xr2,2)*xr4);
b1=6*xr1*xr2*Math.pow(xr3,2)+6*Math.pow(xr1,2)*xr3*xr4;
c1=8*xr1*Math.pow(xr2,2)*xr4+4*Math.pow(xr1,2)*(xr3*xr4+xr2*xr5);
d1=24*Math.pow(xr2,3)*xr3;
xr6=(Math.pow(xr1,3)*(a1+b1+c1+d1)-3*Math.pow(xr1,2)*xr2*xr55)/Math.pow(xr1,6)+xr2*(6/Math.pow(1-pr,4)-xr4)+xr1*(24/Math.pow(1-pr,5)-xr5);

xs1=1/(Math.log(1/(1-ps))*(1-ps));
xs2=xs1*(1/(1-ps)-xs1);
xs3=Math.pow(xs2,2)/xs1+xs1*(1/Math.pow(1-ps,2)-xs2);
xs4=(3*xs1*xs2*xs3-2*Math.pow(xs2,3))/Math.pow(xs1,2)+xs1*(2/Math.pow(1-ps,3)-xs3);
xs5=(-12*xs1*Math.pow(xs2,2)*xs3+3*Math.pow(xs1,2)*Math.pow(xs3,2)+4*Math.pow(xs1,2)*xs2*xs4+6*Math.pow(xs2,4)) /Math.pow(xs1,3)+xs1*(6/Math.pow(1-ps,4)-xs4);

er=xr+pr*qr*xr2/(2*(n + 2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n + 2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3+(1/8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);

let vorder=[];
vorder[0]=er;vorder[1]=crs;
return vorder;
}

//********Нормальные порядковые статистики**************************************

function ordern(n,pr,ps) {

/*
  Computes David-Johnson approximation for mean and covariance between rth
  and sth order statistics from the normal dist. for a sample size n.
  pr=r/(n+1);ps=s/(n+1);
  vorder[0] - mean;
  vorder[1] - covariance between rth and sth order statistics 
  if r=s then vorder[1] - variance
*/
  let qr,qs,xr,xs,dr,ds;
  let xr1,xr2,xr3,xr4,xr5,xr6;
  let xs1,xs2,xs3,xs4,xs5,xs6;
  let z1,z2,z3,z4,z5,z6,z7;
  let er,crs;
  let vorder=[];

 qr=1-pr;qs=1-ps; 
xr=invnormaldistribution(pr);
xs=invnormaldistribution(ps);
dr=s2pi*Math.exp(xr*xr/2.);
ds=s2pi*Math.exp(xs*xs/2.);
xr1=dr;
xr2=xr*Math.pow(dr,2);
xr3=(2*xr*xr+1)*Math.pow(dr,3);
xr4=(6*xr*xr*xr+7*xr)*Math.pow(dr,4);
xr5=(24*Math.pow(xr,4)+46*xr*xr+7)*Math.pow(dr,5);
xr6=(120*Math.pow(xr,5)+326*Math.pow(xr,3)+127*xr)*Math.pow(dr,6);
xs1=ds;
xs2=xs*Math.pow(ds,2);
xs3=(2*xs*xs+1)*Math.pow(ds,3);
xs4=(6*Math.pow(xs,3)+7*xs)*Math.pow(ds,4);
xs5=(24*Math.pow(xs,4)+46*Math.pow(xs,2)+7)*Math.pow(ds,5);
xs6=(120*Math.pow(xs,5)+326*Math.pow(xs,3)+127*xs)*Math.pow(ds,6);

er=xr+pr*qr*xr2/(2*(n+2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n+2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3 + (1 /8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);

vorder[0]=er;
vorder[1]=crs;
return(vorder);
}

//***********************************************************************
function gammad (x,p)


//
//  Purpose:
//
//    GAMMAD computes the Incomplete Gamma Integral
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by B Shea.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    B Shea,
//    Algorithm AS 239:
//    Chi-squared and Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, let X, P, the parameters of the incomplete 
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int IFAULT, error flag.
//    0, no error.
//    1, X < 0 or P <= 0.
//
//    Output, let GAMMAD, the value of the incomplete 
//    Gamma integral.
//
{
  let a;
  let an;
  let arg;
  let b;
  let c;
  let elimit = - 88.0;
  let oflo = 1.0E+37;
  let plimit = 1000.0;
  let pn1;
  let pn2;
  let pn3;
  let pn4;
  let pn5;
  let pn6;
  let rn;
  let tol = 1.0E-14;
   let upper;
  let valuex;
  let xbig = 1.0E+08;

  valuex = 0.0;
//
//  Check the input.
//
  if ( x < 0.0 )
  {
    ifault = 1;
    return valuex;
  }

  if ( p <= 0.0 )
  {
    ifault = 1;
    return valuex;
  }

  ifault = 0;

  if ( x == 0.0 )
  {
    valuex = 0.0;
    return valuex;
  }
//
//  If P is large, use a normal approximation.
//
  if ( plimit < p )
  {
    pn1 = 3.0 * Math.sqrt ( p ) * ( Math.pow ( x / p, 1.0 / 3.0 ) 
    + 1.0 / ( 9.0 * p ) - 1.0 );

    upper = false;
    valuex = 1-normaldistribution(pn1); //alnorm ( pn1, upper );
    return valuex;
  }
//
//  If X is large set valuex = 1.
//
  if ( xbig < x )
  {
    valuex = 1.0;
    return valuex;
  }
//
//  Use Pearson's series Math.expansion.
//  (Note that P is not large enough to force overflow in AMath.logAM).
//  No need to test IFAULT on exit since P > 0.
//
  if ( x <= 1.0 || x < p )
  {
    arg = p * Math.log ( x ) - x - Math.log(gamma(p+1.0));
    c = 1.0;
    valuex = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      valuex = valuex + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + Math.log ( valuex );

    if ( elimit <= arg )
    {
      valuex = Math.exp ( arg );
    }
    else
    {
      valuex = 0.0;
    }
  }
//
//  Use a continued fraction Math.expansion.
//
  else 
  {
    arg = p * Math.log ( x ) - x - Math.log(gamma(p));
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    valuex = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if ( pn6 != 0.0 )
      {
        rn = pn5 / pn6;

        if ( Math.abs ( valuex - rn ) <= r8_min ( tol, tol * rn ) )
        {
          break;
        }
        valuex = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
//
//  Re-scale terms in continued fraction if terms are large.
//
      if ( oflo <= Math.abs ( pn5 ) )
      {
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }

    arg = arg + Math.log ( valuex );

    if ( elimit <= arg )
    {
      valuex = 1.0 - Math.exp ( arg );
    }
    else
    {
      valuex = 1.0;
    }
  }

  return valuex;
}
//****************************************************************************80

function ppchi2 (p,v) {

//****************************************************************************80
//
//  Purpose:
//
//    ppchi2 evaluates the percentage points of the Chi-squared PDF.
//
//  Discussion
//
//    Incorporates the suggested changes in AS R85 (vol.40(1),
//    pages 233-5, 1991) which should eliminate the need for the limited
//    range for P, though these limits have not been removed
//    from the routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Donald Best, DE Roberts.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Donald Best, DE Roberts,
//    Algorithm AS 91:
//    The Percentage Points of the Chi-Squared Distribution,
//    Applied Statistics,
//    Volume 24, Number 3, 1975, pages 385-390.
//
//  Parameters:
//
//    Input, let P,  value of the chi-squared cumulative
//    probability density function.
//    0.000002 <= P <= 0.999998.
//
//    Input, let V, the parameter of the chi-squared probability
//    density function.
//    0 < V.
//
//    Input, let G, the value of Math.log ( Gamma ( V / 2 ) ).
//
//    Output, int *IFAULT, is nonzero if an error occurred.
//    0, no error.
//    1, P is outside the legal range.
//    2, V is not positive.
//    3, an error occurred in GAMMAD.
//    4, the result is probably as accurate as the machine will allow.
//
//    Output, let PPCHI2, the value of the chi-squared random
//    deviate with the property that the probability that a chi-squared random
//    deviate with parameter V is less than or equal to PPCHI2 is P.
//

  let a;
  let aa = 0.6931471806;
  let b;
  let c;
  let c1 = 0.01;
  let c2 = 0.222222;
  let c3 = 0.32;
  let c4 = 0.4;
  let c5 = 1.24;
  let c6 = 2.2;
  let c7 = 4.67;
  let c8 = 6.66;
  let c9 = 6.73;
  let c10 = 13.32;
  let c11 = 60.0;
  let c12 = 70.0;
  let c13 = 84.0;
  let c14 = 105.0;
  let c15 = 120.0;
  let c16 = 127.0;
  let c17 = 140.0;
  let c18 = 175.0;
  let c19 = 210.0;
  let c20 = 252.0;
  let c21 = 264.0;
  let c22 = 294.0;
  let c23 = 346.0;
  let c24 = 420.0;
  let c25 = 462.0;
  let c26 = 606.0;
  let c27 = 672.0;
  let c28 = 707.0;
  let c29 = 735.0;
  let c30 = 889.0;
  let c31 = 932.0;
  let c32 = 966.0;
  let c33 = 1141.0;
  let c34 = 1182.0;
  let c35 = 1278.0;
  let c36 = 1740.0;
  let c37 = 2520.0;
  let c38 = 5040.0;
  let ch;
  let e = 0.5E-06;
  let i;
  let if1;
  let maxit = 20;
  let pmax = 0.999998;
  let pmin = 0.000002;
  let p1;
  let p2;
  let q;
  let s1;
  let s2;
  let s3;
  let s4;
  let s5;
  let s6;
  let t;
  let valuex;
  let x;
  let xx;
  let g;
  g=Math.log(gamma(v/2));
//
//  Test arguments and initialize.
//
  valuex = - 1.0;

  if ( p < pmin || pmax < p )
  {
    ifault = 1;
    return valuex;
  }

  if ( v <= 0.0 )
  {
    ifault = 2;
    return valuex;
  }

  ifault = 0;
  xx = 0.5 * v;
  c = xx - 1.0;
//
//  Starting approximation for small chi-squared
//
  if ( v < - c5 * Math.log ( p ) )
  {
    ch = Math.pow ( p * xx * Math.exp ( g + xx * aa ), 1.0 / xx );

    if ( ch < e )
    {
      valuex = ch;
      return valuex;
    }
  }
//
//  Starting approximation for V less than or equal to 0.32
//
  else if ( v <= c3 )
  {
    ch = c4;
    a = Math.log ( 1.0 - p );

    for ( ; ; )
    {
      q = ch;
      p1 = 1.0 + ch * ( c7 + ch );
      p2 = ch * (c9 + ch * ( c8 + ch ) );

      t = - 0.5 + (c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 
      3.0 * ch ) ) / p2;

      ch = ch - ( 1.0 - Math.exp ( a + g + 0.5 * ch + c * aa ) * p2 / p1) / t;

      if ( Math.abs ( q / ch - 1.0 ) <= c1 )
      {
        break;
      }
    }
  }
  else
  {
//
//  Call to algorithm AS 111 - note that P has been tested above.
//  AS 241 could be used as an alternative.
//
    x = ppnd ( p, ifault );
//
//  Starting approximation using Wilson and Hilferty estimate
//
    p1 = c2 / v;
    ch = v * Math.pow ( x * Math.sqrt ( p1 ) + 1.0 - p1, 3 );
//
//  Starting approximation for P tending to 1.
//
    if ( c6 * v + 6.0 < ch )
    {
      ch = - 2.0 * ( Math.log ( 1.0 - p ) - c * Math.log ( 0.5 * ch ) + g );
    }
  }
//
//  Call to algorithm AS 239 and calculation of seven term
//  Taylor series
//
  for ( i = 1; i <= maxit; i++ )
  {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p - gammad ( p1, xx);
   
    //if ( if1 != 0 )
    //{
    //  ifault = 3;
    //  return valuex;
    //}

    t = p2 * Math.exp ( xx * aa + g + p1 - c * Math.log ( ch ) );
    b = t / ch;
    a = 0.5 * t - b * c;
    s1 = ( c19 + a * ( c17 + a * ( c14 + a * ( c13 + a * ( c12 + 
    c11 * a ))))) / c24;
    s2 = ( c24 + a * ( c29 + a * ( c32 + a * ( c33 + c35 * a )))) / c37;
    s3 = ( c19 + a * ( c25 + a * ( c28 + c31 * a ))) / c37;
    s4 = ( c20 + a * ( c27 + c34 * a) + c * ( c22 + a * ( c30 + c36 * a ))) / c38;
    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a )) / c37;
    s6 = ( c15 + c * ( c23 + c16 * c )) / c38;
    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c * ( s1 - b * 
    ( s2 - b * ( s3 - b * ( s4 - b * ( s5 - b * s6 ))))));

    if ( e < Math.abs ( q / ch - 1.0 ) )
    {
       valuex = ch;
       return valuex;
    }
  }

 ifault = 4;
 valuex = ch;

  return valuex;
}
//****************************************************************************80

function ppnd (p)

//****************************************************************************80
//
//  Purpose:
//
//    PPND produces the normal deviate value corresponding to lower tail area = P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by J Beasley, S Springer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    J Beasley, S Springer,
//    Algorithm AS 111:
//    The Percentage Points of the Normal Distribution,
//    Applied Statistics,
//    Volume 26, Number 1, 1977, pages 118-121.
//
//  Parameters:
//
//    Input, let P, the value of the cumulative probability
//    densitity function.  0 < P < 1.
//
//    Output, integer *IFAULT, error flag.
//    0, no error.
//    1, P <= 0 or P >= 1.  PPND is returned as 0.
//
//    Output, let PPND, the normal deviate value with the property that
//    the probability of a standard normal deviate being less than or
//    equal to PPND is P.
//
{
  let a0 = 2.50662823884;
  let a1 = -18.61500062529;
  let a2 = 41.39119773534;
  let a3 = -25.44106049637;
  let b1 = -8.47351093090;
  let b2 = 23.08336743743;
  let b3 = -21.06224101826;
  let b4 = 3.13082909833;
  let c0 = -2.78718931138;
  let c1 = -2.29796479134;
  let c2 = 4.85014127135;
  let c3 = 2.32121276858;
  let d1 = 3.54388924762;
  let d2 = 1.63706781897;
  let r;
  let split = 0.42;
  let valuex;

  ifault = 0;
//
//  0.08 < P < 0.92
//
  if ( Math.abs ( p - 0.5 ) <= split )
  {
    r = ( p - 0.5 ) * ( p - 0.5 );

    valuex = ( p - 0.5 ) * ( ( ( 
        a3   * r 
      + a2 ) * r 
      + a1 ) * r 
      + a0 ) / ( ( ( ( 
        b4   * r 
      + b3 ) * r 
      + b2 ) * r 
      + b1 ) * r 
      + 1.0 );
  }
//
//  P < 0.08 or P > 0.92,
//  R = min ( P, 1-P )
//
  else if ( 0.0 < p && p < 1.0 )
  {
    if ( 0.5 < p )
    {
      r = Math.sqrt ( - Math.log ( 1.0 - p ) );
    }
    else
    {
      r = Math.sqrt ( - Math.log ( p ) );
    }

    valuex = ( ( ( 
        c3   * r 
      + c2 ) * r 
      + c1 ) * r 
      + c0 ) / ( ( 
        d2   * r 
      + d1 ) * r 
      + 1.0 );

    if ( p < 0.5 )
    {
      valuex = - valuex;
    }
  }
//
//  P <= 0.0 or 1.0 <= P
//
  else
  {
    ifault = 1;
    valuex = 0.0;
  }

  return valuex;
}
//****************************************************************************80

function r8_min (x,y) {

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, let X, Y, the quantities to compare.
//
//    Output, let R8_MIN, the minimum of X and Y.
//

  let valuex;

  if ( y < x )  {
    valuex = y;
  } 
  else  {
    valuex = x;
  }
  return valuex;
}

//**********************************численное интегрирование**********************************
function IntegrateFunction(k,nstep,xl,xu,fpolinom) {

 /*
 k-номер метода и порядок полинома; 3,4-метод Симпсона;5-метода Буля
 nstep-количество точек
 xl-нижний предел
 xu-верхний предел 
 fpolinom-имя интегрируемой функции
*/

  let m,h,i,j,sum,zcoef,sm=[],z;

  switch(k) {
   case "3": //Simpson 1/6
     sm[1] = 1; sm[2] = 4; sm[3] = 1;
     zcoef = 1 / 6;   
     break;
   case "4": //Simpson 1/8
   sm[1] = 1; sm[2] = 3; sm[3] = 3; sm[4] = 1;
     zcoef = 1 / 8;
     break;
   case "5":  //Bool 1/90
     sm[1] = 7; sm[2] = 32; sm[3] = 12; sm[4] = 32; sm[5] = 7;
     zcoef = 1 / 90;
     break;
  }
       sum=0;h=(xu-xl)/nstep;
         for(i=1;i<=nstep;i++) {
            for(j=1;j<=k;j++) {
                z=fpolinom(xl+(i-1)*h+(j-1)*h/(k-1));
                sum=sum+sm[j]*z;
            }
         }
      sum=zcoef*sum*h;
    return(sum);
}
//**нормальная аппроксимация процентных точек нецентрального t-распредления(большие выборки)**
function invnontapp(beta,f,d) {
/*
 Input:
 beta-вероятность
 f-число степеней свободы
 d-параметр нецентральности
 Output: zx-квантиль
*/

 let z,zb,zx,f4x;
 zb=invnormaldistribution(beta);
 f4x=1-1/(4*f);
 z=f4x*f4x-zb*zb/(2.*f);
 if(z<0) {
   z=1;f4x=1;
 }
 zx=(f4x*d+zb*Math.sqrt(z+d*d/(2.*f)))/z;
 return(zx);
}
//*****************************************************************************
function simpl(x,nx,stepx,eps,lim,funx) {

/*
x[nx]–вектор размерности nx,на входе содержащий начальные приближения,
на выходе точки минимума;
nx–число переменных минимизируемой функции;
step–начальный шаг минимизации;
eps – относительная точность выхода;
lim – максимальное число итераций;
simpl возвращает число выполненных итераций;
funx – имя минимизируемой функции
Вызов функции: iter=simpl(x,nx,stepx,eps,lim,funx)
*/

    var x1=[];
    var sum=[];
    var dop=[];
    
    k=0;
    alfasimpl=1;
    betasimpl=0.45;
    gama=2.8;
    istep=0;
    ier=0;
k1=nx+1;
k2=nx+2;
k3=nx+3;
k4=nx+4;
vn=nx;
xnx=1/vn;
step1=stepx/(vn*1.41421356237309)*(Math.pow(vn+1,0.5)+vn-1);
step2=stepx/(vn*1.41421356237309)*(Math.pow(vn+1,0.5)-1);

for(i=2;i<=k1;i++) {
 l=(i-1)*nx;
 l1=l+i-1;
    for(j=1;j<=nx;j++) {
      l2=l+j;
      x1[l2]=x[j]+step2;
    }
    x1[l1]=x[i-1]+step1;
}
for(j=1;j<=nx;j++) {
  x1[j]=x[j];
}

flag="M25";
while(1>0) {
switch(flag) {
case "M25":
for(i=1;i<=k1;i++) {
    l=(i-1)*nx+1;
    for(ik=1;ik<=nx;ik++) {
      dop[ik]=x1[l];
      l++;
    }
      q=funx(dop);
      sum[i]=q;
}
case "M28":
sumh=sum[1];
index=1;
for(i=2;i<=k1;i++) {
  if(sum[i]>sumh) {
   sumh=sum[i];
   index=i;
  }
}
suml=sum[1];
kount=1;
for(i=2;i<=k1;i++) {
  if(suml>sum[i]) {
  suml=sum[i];
  kount=i;
  }
}
istep++;
differ=0;
for(i=1;i<=k1;i++) {
  dif=sum[i]-sum[kount];
  differ+=dif*dif;
}
differ=xnx*(Math.pow(differ,0.5));
if(suml<=eps) {
  if(differ<=eps) {flag="M30"; break;}
  if(differ>eps) {flag="M26"; break;}
}
if((differ/(Math.abs(suml)))<=eps) {flag="M30"; break;}
if((differ/(Math.abs(suml)))>eps) {flag="M26"; break;}
case "M30":
differ=0;
l1=(kount-1)*nx;
for(i=1;i<=k1;i++) {
  if(kount !=i) {
  l=(i-1)*nx;
  for(j=1;j<=nx;j++) {
    l2=l+j;
    l3=l1+j;
    if(x1[l3]<=eps) {
    dif=x1[l2]-x1[l3];
    } else {
    dif=(x1[l2]-x1[l3])/x1[l3];
    }
    differ+=dif*dif;
  }
   }
}
differ=xnx*Math.pow(differ,0.5);
if(differ<=eps) {flag="M23"; break;}
case "M26":
if(istep>=lim) {flag="M38"; break;}
if(k==1) {flag="M17"; break;}
for(j=1;j<=nx;j++) {
  sum2=0;
  for(i=1;i<=k1;i++) {
    l=(i-1)*nx+j;
    sum2+=x1[l];
  }
  l2=(index-1)*nx+j;
  l1=k1*nx+j;
  x1[l1]=(sum2-x1[l2])*xnx;
  l3=l1+nx;
  x1[l3] = (1 + alfasimpl)*x1[l1] - alfasimpl*x1[l2];
}
inx=k2*nx+1;
for(ik=1;ik<=nx;ik++) {
  dop[ik]=x1[inx];
  inx++;
}
    q=funx(dop);
    sum[k3]=q;

if(sum[k3]<suml) {flag="M11"; break;}
sumS=suml;
for(i=1;i<=k1;i++) {
  if(index==i) {
  } else if(sum[i]<=sumS) {
  } else {
  sumS=sum[i];
  }

}
if(sum[k3]>sumS) {flag="M13"; break;}
flag="M14";break;
case "M11":
for(j=1;j<=nx;j++) {
  l=k1*nx+j;
  l1=l+nx;
  l2=l1+nx;
  x1[l2]=(1-gama)*x1[l]+gama*x1[l1];
}
inx=k3*nx+1;
for(ik=1;ik<=nx;ik++) {
  dop[ik]=x1[inx];
  inx++;
}

     q=funx(dop);
      sum[k4]=q;

if(sum[k4]<suml) {flag="M16"; break;}
flag="M14";break;
case "M13":
if(sum[k3]>sumh) {flag="M17"; break;}
k=1;
flag="M14"; break;
case "M17":
for(j=1;j<=nx;j++) {
  l=k3*nx+j;
  l1=(index-1)*nx+j;
  l2=l-nx-nx;
  x1[l]=betasimpl*x1[l1]+(1-betasimpl)*x1[l2];
}
k=0;
inx=k3*nx+1;
for(ik=1;ik<=nx;ik++) {
  dop[ik]=x1[inx];
  inx++;
}

     q=funx(dop);
     sum[k4]=q;

if(sumh>sum[k4]) {flag="M16"; break;}
for(j=1;j<=nx;j++) {
  l1=(kount-1)*nx+j;
  for(i=1;i<=k1;i++) {
    if(i !=kount) {
    l=(i-1)*nx+j;
    }
    x1[l]=0.5*(x1[l]+x1[l1]);
  }
}

for(i=1;i<=k1;i++) {
  if(i !=kount) {
   l=(i-1)*nx+1;
   for(ik=1;ik<=nx;ik++) {
     dop[ik]=x1[l];
      l++;
   }
     q=funx(dop);
     sum[i]=q;
  }
}
flag="M28";
break;
case "M16":
for(j=1;j<=nx;j++) {
  l=(index-1)*nx+j;
  l1=k3*nx+j;
  x1[l]=x1[l1];
  sum[index]=sum[k4];
}
flag="M28";
break;
case "M14":
for(j=1;j<=nx;j++) {
  l=(index-1)*nx+j;
  l1=k2*nx+j;
  x1[l]=x1[l1];
  sum[index]=sum[k3];
}
flag="M28";
break;
case "M38":
ier=1;
case "M23":
for(j=1;j<=nx;j++) {
  l=(kount-1)*nx+j;
  x[j]=x1[l];
}
sum[1]=suml;
lim=istep;

return lim;
} //end switch
} //end while
}
//******************************************************************
