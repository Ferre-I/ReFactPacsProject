
using DelimitedFiles
#using DataStructures
using BenchmarkTools


input0=readdlm(stdin,String)


function topologiacg(input::Array{String})

  natmax=UInt16(10000)
	nmax=UInt16(6000)
	#
	ica=Array{Int64}(zeros(natmax))
	io=Array{Int64}(zeros(natmax))
	ico=Array{Int64}(zeros(natmax))
	inn=Array{Int64}(zeros(natmax))  ##in(natmax)
	ih=Array{Int64}(zeros(natmax))

	rhc=Array{Float64}(zeros(natmax))
	xm=Array{Float64}(zeros(natmax))
	ind1=Array{Int64}(zeros(natmax))
	ind2=Array{Int64}(zeros(natmax))
	inb1=Array{Int64}(zeros(natmax))
	
	##ib
	ibn=Array{Int64}(zeros(natmax))
	ibca=Array{Int64}(zeros(natmax))
	ibco=Array{Int64}(zeros(natmax))
	ibh=Array{Int64}(zeros(natmax))
	ibo=Array{Int64}(zeros(natmax))
	
	
	
	
	xmb=Array{Float64}(zeros(natmax))
	##
	##
	r=Array{Float64}(zeros(natmax,3))
	b=Array{Float64}(zeros(natmax,3))
	ipsb=Array{Float64}(zeros(nmax,nmax))
	##
	rcm=Array{Float32}(zeros(3))
	##
	part=Array{String}(undef,natmax)
	bead=Array{String}(undef,natmax)
	atom=Array{String}(undef,natmax)
	res=Array{String}(undef,natmax)
	cad=Array{String}(undef,natmax)
	cbead=Array{String}(undef,natmax)

	c1=String("")
	c2=String("")
	c3=String("")
	c4=String("")

	natom=Int64(0)
	nres=Int64(0)
	nbead=Int64(0)
	k0=Int64(0)
	kk=Int64(0)
	nres0=Int64(0)

	filas=size(input,1)
	columnas=size(input,2)

	i=Int64(1)
	j=Int64(1)
	k=Int64(1)
	n=Int64(1)
	m=Int64(0)


#=Constantes
=#
	rco=1.4
	rca=2.
	rcb=2.
	rs=1.5
	rn=1.2
	ro=1.2
	rh=0.7

	xmc=.012
	xmo=.016
	xmn=.014
	xms=.032
	xmh=.01
	xmres=.1


	#Filas >=1 and length de fila es igual a 9
	if filas>=0 && columnas==9
		
		c1=input[i,1]
		j=parse(Int64,input[i,2])
		c2=input[i,3]
		c3=input[i,4]
		c4=input[i,5]
		k=parse(Int64,input[i,6])
		r[i,1]=parse(Float64,input[i,7])
		r[i,2]=parse(Float64,input[i,8])
		r[i,3]=parse(Float64,input[i,9])

		atom[i]=c2
		res[i]=c3		
		kk+=1
		cad[i]=c4
		ind2[i]=kk
		k0=k
		ind1[i]=ind2[i]-nres0
		k=ind2[i]	
		
			if atom[i]=="N"  inn[k]=i	end
			if atom[i]=="CA" ica[k]=i	end
			if atom[i]=="C"  ico[k]=i	end
			if atom[i]=="O"  io[k]=i	end
			if atom[i]=="H"  ih[k]=i	end	
	
		
		i+=1
		while i<=filas
			c1=input[i,1]
			j=parse(Int64,input[i,2])
			c2=input[i,3]
			c3=input[i,4]
			c4=input[i,5]
			k=parse(Int64,input[i,6])
			r[i,1]=parse(Float64,input[i,7])
			r[i,2]=parse(Float64,input[i,8])
			r[i,3]=parse(Float64,input[i,9])	
				
				
					
			atom[i]=c2
			res[i]=c3
		
			if c4!=cad[i-1] 
				nres0=ind2[i-1]
			end
			
			if k!=k0
				kk+=1	
			end

			cad[i]=c4
			ind2[i]=kk
			k0=k
			ind1[i]=ind2[i]-nres0
			k=ind2[i]

			if atom[i]=="N"  inn[k]=i	end
			if atom[i]=="CA" ica[k]=i	end
			if atom[i]=="C"  ico[k]=i	end
			if atom[i]=="O"  io[k]=i	end
			if atom[i]=="H"  ih[k]=i 	end
		i+=1
		end	##END WHILE
		n=i
		natom=i-1
###############################################################
	nres=ind2[natom]
	print("$natom $nres\n\n")

		for i=1:natom-1
		xm[i]=xmc
		rhc[i]=rco
			if atom[i]=="CA"
			rhc[i]=rca
			end
			if atom[i]=="CB"
			rhc[i]=rcb
			end
			if atom[i][1:1]=="N" 
			xm[i]=xmn
			rhc[i]=rn
			end
			if atom[i][1:1]=="O"
			xm[i]=xmo
			rhc[i]=ro
			end
			if atom[i][1:1]=="S"
			xm[i]=xms
			rhc[i]=rs
			end
			if atom[i]=="H"
			xm[i]=xmh
			rhc[i]=rh
			end
		end  ##END FOR
		m=nbead+1
		

		for i=1:nres
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			ii=inn[i]
				#print("$ii $m $(cad[ii]) kl\n")
				b[m,1]=r[ii,1]
				b[m,2]=r[ii,2]
				b[m,3]=r[ii,3]
				part[m]="N"
				bead[m]=res[ii]
				cbead[m]=cad[ii]
				inb1[m]=ind1[ii]
				ibn[i]=m
				m=m+1
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			#####
			
			if res[ii]!="PRO"
				ii=ih[i]
				#print("$ii $m $(cad[ii]) hh\n")		
				b[m,1]=r[ii,1]
				b[m,2]=r[ii,2]
				b[m,3]=r[ii,3]
				part[m]="H"
				bead[m]=res[ii]
				cbead[m]=cad[ii]
				inb1[m]=ind1[ii]
				ibh[i]=m
				m=m+1
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			#####
			ii=ica[i]
			b[m,1]=r[ii,1]
			b[m,2]=r[ii,2]
			b[m,3]=r[ii,3]
			part[m]="CA"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			ibca[i]=m
			##$print("IBCA $i $m\n")
			m=m+1
			rcm[1]=rcm[1]+xm[ii]*r[ii,1]
			rcm[2]=rcm[2]+xm[ii]*r[ii,2]
			rcm[3]=rcm[3]+xm[ii]*r[ii,3]
			xmbead=xmbead+xm[ii]
		##
			ii=ico[i]
			b[m,1]=r[ii,1]
			b[m,2]=r[ii,2]
			b[m,3]=r[ii,3]
			part[m]="C"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			ibco[i]=m
			m=m+1
			rcm[1]=rcm[1]+xm[ii]*r[ii,1]
			rcm[2]=rcm[2]+xm[ii]*r[ii,2]
			rcm[3]=rcm[3]+xm[ii]*r[ii,3]
			xmbead=xmbead+xm[ii]
		##
			ii=io[i]
			b[m,1]=r[ii,1]
			b[m,2]=r[ii,2]
			b[m,3]=r[ii,3]
			part[m]="O"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			ibo[i]=m
			rcm[1]=rcm[1]+xm[ii]*r[ii,1]
			rcm[2]=rcm[2]+xm[ii]*r[ii,2]
			rcm[3]=rcm[3]+xm[ii]*r[ii,3]
			xmbead=xmbead+xm[ii]
			xmb[m-2]=xmbead
			m=m+1
		##	
		##Enlaces correspondientes a la cadena lateral		
			if res[ii]=="ARG"
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=1:4
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=5:7
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			end ##END IF
		##
		if res[ii]=="ASN" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=1:4
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			xmb[m]=xmbead
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="ASP"
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=1:4
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		#
		if res[ii][1:2]=="CY"
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=1:2
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="GLN" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:5
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="GLU" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:5
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii][1:2]=="HI" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:2
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=3:4
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=5:6
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S3"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="ILE" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:4
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="LEU" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:4
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="LYS" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=1:3
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
				for j=4:5
					nbead=m-1	
					ii=ica[i]+j
					rcm[1]=rcm[1]+xm[ii]*r[ii,1]
					rcm[2]=rcm[2]+xm[ii]*r[ii,2]
					rcm[3]=rcm[3]+xm[ii]*r[ii,3]
					xmbead=xmbead+xm[ii]
				end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="MET" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:4
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="PHE" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:3
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=4:5
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=6:7
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S3"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="PRO" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:3
				ii=ica[i]-j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="SER" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:2
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="THR" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:3
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="TRP" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:2
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=3:5
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=8:10
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S3"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=6:7
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S4"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="TYR" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:3
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=4:5
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S2"
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			xmb[m]=xmbead
			inb1[m]=ind1[ii]
			m=m+1
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=6:8
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S3"
			xmb[m]=xmbead
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		if res[ii]=="VAL" 
			rcm[1]=0.
			rcm[2]=0.
			rcm[3]=0.
			xmbead=0.
			for j=1:3
				ii=ica[i]+j
				rcm[1]=rcm[1]+xm[ii]*r[ii,1]
				rcm[2]=rcm[2]+xm[ii]*r[ii,2]
				rcm[3]=rcm[3]+xm[ii]*r[ii,3]
				xmbead=xmbead+xm[ii]
			end
			b[m,1]=rcm[1]/xmbead
			b[m,2]=rcm[2]/xmbead
			b[m,3]=rcm[3]/xmbead
			part[m]="S1"
			xmb[m]=xmbead
			bead[m]=res[ii]
			cbead[m]=cad[ii]
			inb1[m]=ind1[ii]
			m=m+1
		end ##END IF
		##
		end	##END FOR##
	###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		nbead=m-1
		for k=0:nres-1
			if k!=0
			i0=ibca[k]
			j0=ibn[k]
			else
			 i0=0
			 j0=0
			end
			
			i1=ibca[k+1]
			j1=ibn[k+1]
			i2=ibca[k+2]
			j2=ibn[k+2]
		  	
			if bead[i1]!="PRO"
				ipsb[i1-2,i1-1]=1
				ipsb[i1-2,i1]=1
				ipsb[i1-2,i1+1]=1
				ipsb[i1-1,i1]=1
				else

				ipsb[i1-1,i1]=1
				ipsb[i1-1,i1+1]=1

				ipsb[i0+1,i1+1]=1
			end
			
			
			ipsb[i1,i1+1]=1
			ipsb[i1,i1+2]=1
			ipsb[i1+1,i1+2]=1
			
			if k<nres -1
			#= Aquí se puso esta estructura para traducir el indexado de Fortran
			a Julia, debido a que FORTRAN maneja indexado con 0, y JUlia maneja
			indexado a partir de 1. Además que Fortran es tolerante a desbordamientos	
		
			  
			  Este problema se dá en la última iteración cuando k==nres-1
			=#
				
				ipsb[i1,i2-1]=1
				ipsb[i1,i2-2]=1
				ipsb[i1,i2]=1
				ipsb[i1+1,i2-2]=1
				ipsb[i1+1,i2-1]=1
				ipsb[i1+1,i2]=1
				ipsb[i1+2,i2-2]=1
				ipsb[i1+2,i2]=1
			end

			imax=j2-1
			if imax>i1+2
				ipsb[j1,i1+3]=1
				ipsb[i1+1,i1+3]=1
			end
			if k==nres-1 imax=nbead end
			for i=i1+3:imax
			ipsb[i1,i]=1
				for j=i+1:imax
				ipsb[i,j]=1
				end
			end
			if part[i1+3]=="S1"
				if bead[k]!="PRO" 
				ipsb[i1-2,i1+3]=1
				else
				ipsb[i1-1,i1+3]=1
				end
				ipsb[i1+1,i1+3]=1
			end
		end
	###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	for i=1:natom-1
		if atom[i]=="SG"
			for j=i+1:natom
				if atom[j]=="SG"
					rij1=r[i,1]-r[j,1]
					rij2=r[i,2]-r[j,2]
					rij3=r[i,3]-r[j,3]
					rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
						if rij<2.5
							ii=ind2[i]
							jj=ind2[j]
							i1=ibca[ii]
							i2=ibca[jj]
							ipsb[i1+3,i2+3]=1
							ipsb[i1+3,i2]=1
							ipsb[i1,i2+3]=1
						end
				end
			end
		end
	end
	#####
	c1="ATOM"
	
	print("natom0,natom\t0\t$natom\t$nbead\n")
	
	if nbead > nmax
		print("Masa de partículas! $nbead\n")
		print("El maxim possible es $nmax")
	exit(0)
	end
	
	print("0  $natom\n")
	
	wre= open("topcgJulia.dat","w")
		for i=1:nbead-1
			for j=i+1:nbead
				rij1=b[i,1]-b[j,1]
				rij2=b[i,2]-b[j,2]
				rij3=b[i,3]-b[j,3]
				rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
				
				if cbead[i]==cbead[j]
				  if ipsb[i,j]==1 write(wre,"    $(lpad(i,8))    $(lpad(j,8))   $rij\n") end
				end			    #           #         
			end
		end
	close(wre)


		c1="ATOM"
		l=1
		wre=open("structurecgJulia.pdb","w")
		for i=1:nbead
			write(wre, "$c1   $(lpad(i,4))  $(rpad(part[i],2))  $(bead[i]) $(cbead[i]) $(lpad(inb1[i],3))    $(SubString(lpad(b[i,1],8),1,8)) $(SubString(lpad(b[i,2],8),1,8)) $(SubString(lpad(b[i,3],8),1,8))\n")
		end
		close(wre)


	else
	   #println("\n##Error el archivo de entrada no es el formato de salida de gentop\n")
	end
	##FIN IF	
end


@time topologiacg([""])
@time topologiacg([""])
@time topologiacg([""])


topologiacg(input0)
@time topologiacg(input0)

