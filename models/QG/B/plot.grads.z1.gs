function test(args)

'open qgmodelT21.ctl'

'set x 1 64'
'set y 1 32'
'set t 1 1200'

'set fwrite -sq psidiff.z1.dat' 
'set gxout fwrite'
'set z 1 3'
'define fdif=psi(t+1)-psi(t+0)'
'set z 1'
'define ldif=psi(t+1)-psi(t+0)'

ww=0
yy=1
while(yy<33)
xx=1
while(xx<65)
zz=1
while(zz<4)

'set t 1 1200'

'set x ' xx
'set y ' yy
'set z ' zz
'define ddif=fdif'

'set x 1 64'
'set y 1 32'
'define pdif=ddif*ldif'
'set t 1'
ww=ww+1
'd ave(pdif,t=1,t=1200)'
*say "X:" xx "  Y:" yy "  Z: " zz 
*'q define'
*say result 
*say "     "
say "ON WRITEOUT NUMBER: " ww "X:" xx "  Y:" yy "  Z: " zz  
zz=zz+1
endwhile
xx=xx+1   
endwhile
yy=yy+1
endwhile

say "TOTAL WRITES: " ww
'quit'
return
