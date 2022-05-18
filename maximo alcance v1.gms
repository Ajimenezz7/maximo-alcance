**************************************
* Antonio José Jiménez-Zarza Carriazo*
**************************************

$set     n 100

Set      j /0*%n%/;

Sets
         jlast(j),
         jnotlast(j);
         jlast(j)$(ord(j)=card(j))=yes;
         jnotlast(j)=not jlast(j);

Scalar


         h0      'Altura inicial en [m]'                            /2773/,
         h1      'Altura final en [m]'                              /0/,

         x0      'Posicion inicial en [m]'                          /0/,

         rho0    'Densidad en [Kg/m^3]'                             /1.225/,
         g0      'Gravedad en [m/s^2]'                              /9.81/,
         K       'Coeficiente de resistencia inducida [Adim]'       /0.0573/,
         Cd0     'Coeficiente de Resistencia parásita[Adim]'        /0.0209/,
         jota    'Coeficiente de resistenica [Adim]'                /0.0464/,
         S       'Superficie en[m^2]'                               /8/,
         RT      'Radio de la tierra en [m]'                        /6378000/,
         Lambda  'Gradiente de Temperatura [K/m]'                   /0.0065/,
         R       'Constante Gases Ideales'                          /287.57/,
         T0      'Temperatura Nivel del Mar en [K]'                 /288.15/,



         Clmin   'Coeficiente de sustentacion minimo'               /0/,
         Clmax   'Coeficiente de sustentacion maximo'               /0.7347/,

         Vmin    'Velocidad Minima en [m/s]'                        /0/,
         Vmax    'Velocidad Maxima en [m/s]'                        /100/,

         MinGamma_d 'Minimo Angulo Descenso en [rad]'               /0/,
         MaxGamma_d 'Maximo Angulo Descenso en [rad]'               /0.1/

         m 'masa del planeador [kg]'                                /245/

         v10 'velocidad del viento a 10 metros [m/s]'               /30/
         a 'exponente de Hellmann [Adim]'                           /0.14/

         AnguloViento   'Ángulo Viento en [rad] (Viento Descendente)'  /8/;


Parameter Theta 'Conversión de radianes a grados del ángulo del viento';

         Theta = AnguloViento*(pi/180);


Variable

         Cl(j)         'Coeficiente de Sustentación [Adim]',
         h(j)          'Altura [m]',
         V(j)          'Velocidad [m/s]',
         x(j)          'Posicion [m]',
         step          'Escalon tiempo (UNIFORME) [s]',
         Gamma_d(j)    'Angunlo de descenso [rad]',
         Xmax          'Posición [m]',
         g(j)          'Gravedad [m/s^2]',
         rho(j)        'Densidad [Kg/m^3]'
         Vw(j)         'Velocidad del viento [m/s]'

;

Equations

         Rest0,
         Rest1,
         Rest2,
         Rest3,
         Rest4,
         Rest5,
         Rest6,
         Fobj;


* ECUACIÓN DEL VIENTO:

Rest0(j)..       Vw(j) =e=(-1*10**(-12)*(h(j))**4)
                                + (7*10**(-9)*(h(j))**3)
                                         + (-2*10**(-5)*(h(j))**2)
                                                  + (0.017*h(j))+9.2226;


* ECUACIÓN DE LA DENSIDAD:

Rest1(j)..                        rho(j) =e= rho0*(1-(lambda*h(j))/
                                             (T0))**((g0)/(R*lambda)-1);

* ECUACIÓN DE LA GRAVEDAD:

Rest2(j)..                        g(j) =e= g0*((RT*RT)/(RT*RT+h(j)*h(j)));

*ECUACIÓN DINAMICA DE LA SUSTENTACION

Rest3(j)..                        0.5*rho(j)*S*V(j)*V(j)*Cl(j) =e=
                                         m*g(j)*cos(Gamma_d(j));

*ECUACIÓN DINAMICA DE LA RESISTENCIA

Rest4(j)..                        0.5*rho(j)*S*V(j)*V(j)*(Cd0+K*Cl(j)*Cl(j)) =e=
                                         m*g(j)*sin(Gamma_d(j));

*ECUACIÓN CINEMATICA DE LA ALTURA


Rest5(j)$(jnotlast(j))..          h(j+1)-h(j)=e=0.5*(-V(j)*sin(Gamma_d(j))
                                         -V(j+1)*sin(Gamma_d(j+1))-
                                                 Vw(j)*sin(Theta)-
                                                 Vw(j+1)*sin(Theta))*step;

*ECUACION CINEMATICA DE LA VELOCIDAD


Rest6(j)$(jnotlast(j))..          x(j+1)-x(j)=e=0.5*(V(j)*cos(Gamma_d(j))
                                         +V(j+1)*cos(Gamma_d(j+1))-
                                                 Vw(j)*cos(Theta)-
                                                 Vw(j+1)*cos(Theta))*step;

*FUNCION OBJETIVO

Fobj..                            xmax =e= x('%n%');

* BOUNDARY CONDITIONS

         Cl.lo(j) = Clmin;
         Cl.up(j) = Clmax;

         V.lo(j) = vmin;
         V.up(j) = vmax;

         g.lo(j) = 0;
         g.up(j) = g0;

         Gamma_d.lo(j) = MinGamma_d;
         Gamma_d.up(j) = MaxGamma_d;

         rho.lo(j) = 0;
         rho.up(j) = rho0;

* VALORES INICIALES Y FINALES

         h.fx('0')   = h0;
         h.fx('%n%') = h1;
         x.fx('0')   = x0;
         V.fx('0') =  32.9481;

display m;

Model MaximumGliderRange /all/;

Option NLP=ipopt;

Solve MaximumGliderRange using NLP maximizing xmax;

execute_unload "results.gdx" v

execute 'gdxxrw.exe results.gdx o=results.xls var=v.L'

