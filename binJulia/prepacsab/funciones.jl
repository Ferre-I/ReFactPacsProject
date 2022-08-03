function potencial(natom)

	rsolv2=rsolv*rsolv
	rbox2=0.5*rbox
	nres=ind2[natom]
	dwatd=dwat/sqrt(3.0)
	

	 for i=1:natom		#for principal 1
		r1x=r[i,1]+dwat
		r2x=r[i,1]-dwat
		r3x=r[i,1]
		r4x=r[i,1]
		r5x=r[i,1]
		r6x=r[i,1]
		r1y=r[i,2]
		r2y=r[i,2]
		r3y=r[i,2]+dwat
		r4y=r[i,2]-dwat
		r5y=r[i,2]
		r6y=r[i,2]
		r1z=r[i,3]
		r2z=r[i,3]
		r3z=r[i,3]
		r4z=r[i,3]
		r5z=r[i,3]+dwat
		r6z=r[i,3]-dwat
		rv1x=r[i,1]+dwatd
		rv2x=r[i,1]+dwatd
		rv3x=r[i,1]+dwatd
		rv4x=r[i,1]+dwatd
		rv5x=r[i,1]-dwatd
		rv6x=r[i,1]-dwatd
		rv7x=r[i,1]-dwatd
		rv8x=r[i,1]-dwatd
		rv1y=r[i,2]+dwatd
		rv2y=r[i,2]+dwatd
		rv3y=r[i,2]-dwatd
		rv4y=r[i,2]-dwatd
		rv5y=r[i,2]+dwatd
		rv6y=r[i,2]+dwatd
		rv7y=r[i,2]-dwatd
		rv8y=r[i,2]-dwatd
		rv1z=r[i,3]+dwatd
		rv2z=r[i,3]-dwatd
		rv3z=r[i,3]+dwatd
		rv4z=r[i,3]-dwatd
		rv5z=r[i,3]+dwatd
		rv6z=r[i,3]-dwatd
		rv7z=r[i,3]+dwatd
		rv8z=r[i,3]-dwatd
		icont[i]=0
			 for j=1:natom    #for anidado 1
				if imol[j]==imol[i]
					if j!=i
						rij1=r1x-r[j,1]
						  if rij1>rbox2 
						  	rij1=rij1-rbox
						  else rij1<-rbox2
						  rij1=rij1+rbox
						  end
						rij2=r1y-r[j,2]
						  if rij2 > rbox2 
						  rij2=rij2-rbox
						  else rij2 < -rbox2 
						  rij2=rij2+rbox
						  end
						rij3=r1z-r[j,3]
						  if rij3 > rbox2 
						  rij3=rij3-rbox
						  else rij3<-rbox2 
						  rij3=rij3+rbox
						  end
						rmod2=rij1*rij1+rij2*rij2+rij3*rij3
						if rmod2<rsolv2 
							icont[i]=icont[i]+1
							break
						end
					end
				end
			end      #for anidado 1
		 for j=1:natom	  #for anidado 2
			  if imol[j]==imol[i]
					if j!=i
						  rij1=r2x-r[j,1]
						  if rij1 > rbox2 
						  	rij1=rij1-rbox
						  else rij1<-rbox2 
						  	rij1=rij1+rbox
						  end
						  rij2=r2y-r[j,2]
						  if rij2 > rbox2 
						  	rij2=rij2-rbox
						  else  rij2<-rbox2 
						  	rij2=rij2+rbox
						  end
						  rij3=r2z-r[j,3]
						  if rij3>rbox2 
						  	rij3=rij3-rbox
						  else rij3<-rbox2 
						  	rij3=rij3+rbox
						  end
						  rmod2=rij1*rij1+rij2*rij2+rij3*rij3
						if rmod2<rsolv2 
							icont[i]=icont[i]+1
							break
						end
					end
			end
		end	#for anidado 2
		 for j=1:natom		#for anidado 3
			   if imol[j]==imol[i]
					if j!=i
						  rij1=r3x-r[j,1]
						  if rij1>rbox2 
						  	rij1=rij1-rbox
						  else rij1<-rbox2 
						       rij1=rij1+rbox
						  end
						  rij2=r3y-r[j,2]
						  if rij2 > rbox2 
						  	rij2=rij2-rbox
						  else rij2<-rbox2 
						  	rij2=rij2+rbox
						  end
						  rij3=r3z-r[j,3]
						  if rij3>rbox2 
						  	rij3=rij3-rbox
						  else rij3<-rbox2 
						  	rij3=rij3+rbox
						  end
						  rmod2=rij1*rij1+rij2*rij2+rij3*rij3
						  if rmod2<rsolv2
							icont[i]=icont[i]+1
							break
						  end
					end
			end
		end		#for anidado 3
		
		 for j=1:natom	#for anidado 4
			if imol[j]==imol[i]
				if j!=i
					rij1=r4x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=r4y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=r4z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					if rmod2<rsolv2
					icont[i]=icont[i]+1
					break
					end
				end
			end		
		end #for anidado 4
	 for j=1:natom	    #for anidado 5
			if imol[j]==imol[i]
				if j!=i
					rij1=r5x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=r5y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=r5z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2 
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					  if rmod2<rsolv2 
						icont[i]=icont[i]+1
						break
					  end
				end
			end
	end #for anidado 5

	 for j=1:natom #for anidado 6
		if imol[j]==imol[i]
			if j!=i
				rij1=r6x-r[j,1]
				  if rij1>rbox2
				  rij1=rij1-rbox
				  else rij1<-rbox2 
				  rij1=rij1+rbox
				  end
				rij2=r6y-r[j,2]
				  if rij2>rbox2
				  rij2=rij2-rbox
				  else rij2<-rbox2 
				  rij2=rij2+rbox
				  end
				rij3=r6z-r[j,3]
				  if rij3>rbox2
				  rij3=rij3-rbox
				  else rij3<-rbox2 
				  rij3=rij3+rbox
				  end
				rmod2=rij1*rij1+rij2*rij2+rij3*rij3
				if rmod2<rsolv2 
				icont[i]=icont[i]+1
				break
				end
			end
		end
	end #for anidado 6
#!C vertexs
	 for j=1:natom  #for anidado 7
		if imol[j]==imol[i]
			if j!=i
					rij1=rv1x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=rv1y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=rv1z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2 
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					if rmod2<rsolv2 
					icont[i]=icont[i]+1
					break
					end
			end
		end
	end #for anidado 7
	
	 for j=1:natom	#for anidado 8
		if imol[j]==imol[i]
			if j!=i
				rij1=rv2x-r[j,1]
				  if rij1>rbox2
				  rij1=rij1-rbox
				  else rij1<-rbox2 
				  rij1=rij1+rbox
				  end
				rij2=rv2y-r[j,2]
				  if rij2>rbox2
				  rij2=rij2-rbox
				  else rij2<-rbox2 
				  rij2=rij2+rbox
				  end
				rij3=rv2z-r[j,3]
				  if rij3>rbox2
				  rij3=rij3-rbox
				  else rij3<-rbox2 
				  rij3=rij3+rbox
				  end
				rmod2=rij1*rij1+rij2*rij2+rij3*rij3
				if rmod2<rsolv2 
				icont[i]=icont[i]+1
				break
				end
			end
		end
	end #for anidado 8
	
	 for j=1:natom	#for anidado 9
		if imol[j]==imol[i]
			if j!=i
				rij1=rv3x-r[j,1]
				  if rij1>rbox2
				  rij1=rij1-rbox
				  else rij1<-rbox2 
				  rij1=rij1+rbox
				  end
				rij2=rv3y-r[j,2]
				  if rij2>rbox2
				  rij2=rij2-rbox
				  else rij2<-rbox2 
				  rij2=rij2+rbox
				  end
				rij3=rv3z-r[j,3]
				  if rij3>rbox2
				  rij3=rij3-rbox
				  else rij3<-rbox2 
				  rij3=rij3+rbox
				  end
				rmod2=rij1*rij1+rij2*rij2+rij3*rij3
				if rmod2<rsolv2 
				icont[i]=icont[i]+1
				break
				end
			end
		end
	end #for anidado 9
	 for j=1:natom	#for anidado 10
		if imol[j]==imol[i]
			if j!=i
					rij1=rv4x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=rv4y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=rv4z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2 
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					if rmod2<rsolv2 
					icont[i]=icont[i]+1
					break
					end
			end
		end
	end #for anidado 10
	 for j=1:natom	#for anidado 11
		if imol[j]==imol[i]
			if j!=i
					rij1=rv5x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=rv5y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=rv5z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2 
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					if rmod2<rsolv2 
					icont[i]=icont[i]+1
					break
					end
			end
		end
	end #for anidado 11
	
	 for j=1:natom	#for anidado 12
			if imol[j]==imol[i]
				if j!=i
					rij1=rv6x-r[j,1]
					  if rij1>rbox2
					  rij1=rij1-rbox
					  else rij1<-rbox2 
					  rij1=rij1+rbox
					  end
					rij2=rv6y-r[j,2]
					  if rij2>rbox2
					  rij2=rij2-rbox
					  else rij2<-rbox2 
					  rij2=rij2+rbox
					  end
					rij3=rv6z-r[j,3]
					  if rij3>rbox2
					  rij3=rij3-rbox
					  else rij3<-rbox2 
					  rij3=rij3+rbox
					  end
					rmod2=rij1*rij1+rij2*rij2+rij3*rij3
					if rmod2<rsolv2 
					icont[i]=icont[i]+1
					break
					end
				end
			end
	end #for anidado 12
	 for j=1:natom		#for anidado 13
		if imol[j]==imol[i]
			if j!=i
				rij1=rv7x-r[j,1]
				  if rij1>rbox2
				  rij1=rij1-rbox
				  else rij1<-rbox2 
				  rij1=rij1+rbox
				  end
				rij2=rv7y-r[j,2]
				  if rij2>rbox2
				  rij2=rij2-rbox
				  else rij2<-rbox2 
				  rij2=rij2+rbox
				  end
				rij3=rv7z-r[j,3]
				  if rij3>rbox2
				  rij3=rij3-rbox
				  else rij3<-rbox2 
				  rij3=rij3+rbox
				  end
				rmod2=rij1*rij1+rij2*rij2+rij3*rij3
				if rmod2<rsolv2 
				icont[i]=icont[i]+1
				break
				end
			end
		end
	end #for anidado 13
	 for j=1:natom		#for anidado 14
		if imol[j]==imol[i]
			if j!=i
				rij1=rv8x-r[j,1]
				  if rij1>rbox2
				  rij1=rij1-rbox
				  else rij1<-rbox2 
				  rij1=rij1+rbox
				  end
				rij2=rv8y-r[j,2]
				  if rij2>rbox2
				  rij2=rij2-rbox
				  else rij2<-rbox2 
				  rij2=rij2+rbox
				  end
				rij3=rv8z-r[j,3]
				  if rij3>rbox2
				  rij3=rij3-rbox
				  else rij3<-rbox2 
				  rij3=rij3+rbox
				  end
				rmod2=rij1*rij1+rij2*rij2+rij3*rij3
				if rmod2<rsolv2 
				icont[i]=icont[i]+1
				break
				end
			end
		end
	end #for anidado 14

	end #for principal 1
	
	 for i=1:natom   #for principal 2
	fcont[i]=1.0/(1.0+exp((icont[i]-asolv)/bsolv))
	end #for
	
	 for i=1:natom-1 #for principal 3
	ii=ind2[i]
		 for j=i+1:natom  #for anidado 1
		jj=ind2[j]
			if icov[i,j]==0
				if istruct[i,j]==0
					rvdwij=rvdw[i]+rvdw[j]
					sto=(2.0/(nat[i]^0.33+nat[j]^0.33))^6
					potvdw=sqrt(evdw[i]*evdw[j])*sto*(sto-2.0)
					potlk=-0.09/xlamb*(gfree[i]*vol[j]+gfree[j]*vol[i])/(rvdwij^2*exp((rvdwij/xlamb)^2))
					eij=fvdw*potvdw+fsolv*potlk*fcont[i]*fcont[j]+eps*qq[i]*qq[j]/rvdwij
					nstep[i,j]=2
					   rstep[i,j,1]=0.9*rvdwij
					   rstep[i,j,2]=1.1*rvdwij
					 if eij<0.0 
					   estep[i,j,1]=3.0*eij
					   estep[i,j,2]=-eij
					  else
					   estep[i,j,1]=-eij
					   estep[i,j,2]=-eij
					 end
				end
			end
		end #for anidado 1
	end #for principal 3
	
	return
	end
#!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function enchufa(natom,dcut)
	
	dcut2=dcut*dcut
	for  i=1:natom-1
#!C	for  j=i+1:natom
		for  l=1:ipot[i]
		j=npot[i,l]
	#!C	if[istruct[i,j].eq.0]then
	#!C	if[icov[i,j].eq.0]then
			if nat[i]>1&&nat[j]>1
		#!C	  inter[i,j]=0
		#!C mira si hi ha definits potencials d'estructura
			  rij1=dbox[i,j,1]
			  rij2=dbox[i,j,2]
			  rij3=dbox[i,j,3]
			  rmod2=rij1*rij1+rij2*rij2+rij3*rij3
			  if rmod2<dcut2 inter[i,j]=1
			  end
			end
	#!C	endif
	#!C	endif
		end
	end

	nres=ind2[natom]

	for  i=1:nres-1
		n1=in[i]
		n2=in[i+1]
		inter[n1,n2]=0
	end
	for  i=2:nres
		n1=ico[i-1]
		n2=ico[i]
		inter[n1,n2]=0
	end

end


##----------------------------------------------------------------------------

function dmdshake(natom,nbound,mem1,mem2)
	
	ierr=0

#!C particules NO enllaçades
	for i=1:natom-1
#!C	if[i.eq.mem1.or.i.eq.mem2]cycle
		for l=1:ishk[i]
		j=nshk[i,l]
	#!C	if[j.eq.mem1.or.j.eq.mem2]cycle
		rij1=dbox[j,i,1]
		rij2=dbox[j,i,2]
		rij3=dbox[j,i,3]
		rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
		vij1=v[i,1]-v[j,1]
		vij2=v[i,2]-v[j,2]
		vij3=v[i,3]-v[j,3]
		prod=rij1*vij1+rij2*vij2+rij3*vij3
		dmin=rhc[i]+rhc[j]
	#!C xoc frontal entre particules no enllaçades
			if rij < dmin && prod < 0
				chgmom(i,j,rij1,rij2,rij3)
				ierr=ierr+1
			end
		end     #end for
	end

#!C particules enllaçades
	for k=1:nbound
		i=ibound[k,1]
		j=ibound[k,2]
	#!C	if[i.eq.mem1.or.i.eq.mem2]cycle
	#!C	if[j.eq.mem1.or.j.eq.mem2]cycle
		rbmin=rbound[k]*(1.0-sigma)
		rbmax=rbound[k]*(1.0+sigma)
		rbmin2=rbmin*rbmin
		rbmax2=rbmax*rbmax
		rij1=dbox[j,i,1]
		rij2=dbox[j,i,2]
		rij3=dbox[j,i,3]
		rmod2=rij1*rij1+rij2*rij2+rij3*rij3
		vij1=v[i,1]-v[j,1]
		vij2=v[i,2]-v[j,2]
		vij3=v[i,3]-v[j,3]
		prod=rij1*vij1+rij2*vij2+rij3*vij3
		if rmod2 > rbmax2 && prod > 0
		ierr=ierr+1
		chgmom(i,j,rij1,rij2,rij3)
		else
		  if rmod2 < rbmin2 && prod < 0 
		  ierr=ierr+1
		  chgmom(i,j,rij1,rij2,rij3)
		  end
		end
	end

end


#----------------------------------------------------------------------



function chgmom(mem1,mem2,rij1,rij2,rij3)	
	
	
	vdmod=0
	vdmod=vdmod+(v[mem2,1]-v[mem1,1])*rij1
	vdmod=vdmod+(v[mem2,2]-v[mem1,2])*rij2
	vdmod=vdmod+(v[mem2,3]-v[mem1,3])*rij3
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	vdmod=vdmod/rmod2
	xsum=0.5*(1.0/xm[mem1]+1.0/xm[mem2])
#C modul del moment transferit en la colisio

	dp=vdmod/xsum
        v[mem1,1]=v[mem1,1]+dp/xm[mem1]*rij1
        v[mem1,2]=v[mem1,2]+dp/xm[mem1]*rij2
        v[mem1,3]=v[mem1,3]+dp/xm[mem1]*rij3
        v[mem2,1]=v[mem2,1]-dp/xm[mem2]*rij1
        v[mem2,2]=v[mem2,2]-dp/xm[mem2]*rij2
        v[mem2,3]=v[mem2,3]-dp/xm[mem2]*rij3
        
        
	end

##	

	function chgmomene(mem1,mem2,rij1,rij2,rij3,dpot,ich)
	a=1.0e-10
	vdmod=0
	vdmod=vdmod+[v[mem2,1]-v[mem1,1]]*rij1
	vdmod=vdmod+[v[mem2,2]-v[mem1,2]]*rij2
	vdmod=vdmod+[v[mem2,3]-v[mem1,3]]*rij3
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	
#!C projeccio del moment en l'eix que uneix les dues particules
	vdmod=vdmod/rmod2
	xsum=0.5*(1.0/xm[mem1]+1.0/xm[mem2])
#!C modul del moment transferit/distancia en un xoc elastic
	dp=vdmod/xsum
	sto=dp*dp/4.0-dpot/(rmod2*xsum*a*a)
		if sto>0
	#!C sempre es la resta dels dos valors absoluts
			if vdmod > 0  
			dp=dp/2.0-sqrt(sto)
			else
			dp=dp/2.0+sqrt(sto)
			end
		ich=1
		else
	#!C no traspassa la barrera
		ich=0
        end
	v[mem1,1]=v[mem1,1]+dp/xm[mem1]*rij1
        v[mem1,2]=v[mem1,2]+dp/xm[mem1]*rij2
        v[mem1,3]=v[mem1,3]+dp/xm[mem1]*rij3
        v[mem2,1]=v[mem2,1]-dp/xm[mem2]*rij1
        v[mem2,2]=v[mem2,2]-dp/xm[mem2]*rij2
        v[mem2,3]=v[mem2,3]-dp/xm[mem2]*rij3
	end
#!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


	
function creapouhb(n1,n2,rmin,r0,rmax,hb)


#=
Agrega información de una pareja de átomos n y n-1
dentro de una proteína. Con base en el índice n.


Variables globales
inter,istruct,rstep,estep

        ehb es constante igual a 3.

	estep[n1,n2,1]=-1.5*ehb
	estep[n1,n2,2]=1.5*ehb
=#
	inter[n1,n2]=1
	istruct[n1,n2]=1
	nstep[n1,n2]=2
	rstep[n1,n2,1]=rmin
	rstep[n1,n2,2]=rmax
	estep[n1,n2,1]=-4.5
	estep[n1,n2,2]=4.5
end

#-------------------------------------------------------------------------------------------

function dbox(n1,n2,k)

	
#=Variables globales
	rbox,rbox2,r[n1,k]
	
	n:índice de atomo
	k: índice del eje en {x,y,z}
	
	r12 es la distancia entre un átomo n y n-1
	en una proteína
	
	rbox ¡Investigar!
	
	Al final dbox<-r12

Variable local 
r12

=#
	rbox2=0.5*rbox
	r12=r[n2,k]-r[n1,k]
	if r12 > rbox2 then
		r12=r12-rbox
	elseif r12 < -rbox2
		r12=r12+rbox
	end
	dbox=r12

end

#--------------------------------------------------------------------------

  function rnd_gauss(fi, xm, t)
  
  	#Utiliza número pseudoaleatorios
  	##Revisar la generación de la semilla en el código original
  
  	#=
  	
  	fi es una variable generada en la función principal
  	en un principio es un número aleatorio que se guarda en v[i,j]
  	
  	=#
  
  
	R=8.314
	
        std_dev = sqrt((t*R) / xm)
        
        rnd1= rand(Float64)
        rnd2= rand(Float64)
        
        fi = std_dev * sqrt(-2.0 * log( rnd1 ) ) * cos( 2.0 * π * rnd2 )
      end


#----------------------------------------------------------------------------



















