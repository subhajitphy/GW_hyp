**The source code for plotting every element is ```gw_waveform_res.py```**

All usage has been shown in the enclosed jupyter notebook.

* M in total mass in solar mass.
* q is mass ratio.
* et0 is is the initial eccentricity.
* b is the impact parameter in the units of M. 
* ti and tf are initial and final epoch in years.
* t_step the the number of points to be considered while sampling total time of evolution during encounter. 
* inc is the inclination angle.
* distance is the distance to the source (BBH) in Gpc. 
* psrra is the PSR RA and Dec in radian. similar for GW source. 


The $h_{+,\times}$, pre and post-fir residuals can be obtined with
```z=waveform(M,q,et0,b,ti,tf,t_step,inc,distance,psrra, psrdec,gwra,gwdec,rrmethod='None') ```

* z.hp, z.hx, z.response, z.prefitres, z.postfitres returns the hp, hx array,   
reponse of GW (F_+ h_+ +F_x h_x), prefit residuals and post fit residuals respectively. 

* If rrmethod='dudt' the solver used usual LAL implement : (u,e_t,n) approach. BY DEFAULT it is set to be in accordance with PTA approach: (l,e_t,n).

```z1=waveform(M,q,et0,b,ti,tf,t_step,inc,distance,estimatepeak='True')```

z1.hp, z1.hx, z1.x, z1.y, z1.peak returns the $h_+$, $h_\times$ array, [orbit] and peak frequency (in Hz) respectively. 

