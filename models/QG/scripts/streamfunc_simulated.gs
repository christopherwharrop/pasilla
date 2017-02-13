function test(args)

resol=subwrd(args,1)
lev=subwrd(args,2)
expid=subwrd(args,3)

offset=0
* light grey (reserve 81-89)
'set rgb 81 180 180 180'
'set rgb 82 160 160 160'

'open qgmodelT'resol'.ctl'

'q file'
line=sublin(result,5)
ntime=subwrd(line,12)

'set lon -180 180'
'set lat 0 90'
'set lev 'lev
'set vpage 0 11 0 8.5'
'set mproj nps'

'set grads off'
'set font 0'
'set gxout contour'

'set xlopts 1 3 0.18'
'set ylopts 1 3 0.18'

'set clevs 9e32'
'set mpdraw off'
'set grid off'
'd psi'
'run basemap.gs L 81 81 M'

'set ccolor 1'
'set cthick 7'
'set clopts 1 1.2 0.1'
*'set clab off'
'set mpdraw off'
'set grid on 5 82 10'
'set poli off'
'set xaxis 0 360 30'
'set ylevs 15 30 45 60 75'
'd 1e-7*ave(psi,t=1,t='ntime')'
'set string 1 tc'
'set strsiz 0.16 0.18'
'draw string 5.5 8. Mean T'resol'-'expid' streamfunction at 'lev

'enable print st.gx'
'print'
'disable print'
'!gxeps -c st.gx'
'!epstopdf st.eps'
'!mv st.pdf T'resol'-'expid'.mean'lev'.pdf'
'!mv st.eps T'resol'-'expid'.mean'lev'.eps'
'quit'
return
