# I check the velocity, it's about 2.5km/s 0.3Hz-1Hz

using Geodesics
using Dates
using HDF5
using Mmap
using FFTW
using Statistics
using SAC

const LOCATION_FILE="/Volumes/yy2/DATA/locations.txt"
const FFT_DIR="/Volumes/yy2/DATA/FFT/2012/"
const CCF_DIR="/Users/zma/v2.0_4Xperformace/testing.dir/2012_forSSA_C3correlate_BHENZ_daily/"

const d1 = Date(2012,1,1)
const d2 = Date(2012,12,30)
const nseg = 47
const nt_wrong = 9112
const nt = 9113
const nt_original=18225
const downsamp_freq=5.0

const n1 = Int32(200*downsamp_freq)    #this is about 2*200km/(2km/s)*20
const nlen = Int32(1200*downsamp_freq+1) #window length for c3, easier if it's odd so time 0 is centered
const n2 = n1+nlen-1
const n3 = div(nlen,2)+1

const comp_s = ["BHE" "BHN" "BHZ"]
const comp_r = ["BHZ"]
const comp_a = ["BHE" "BHN" "BHZ"]

#const stalist = ["BLOW" "GNW" "STOR" "HOOD" "TOLT"]
const stalist = ["BLOW" "DOSE" "FISH" "GNW" "HOOD" "LEBA" "LON" "RATT" "SP2" "STOR"]

# for fft1*fft2, ifft, fft, c3, load data
const timeall = [0.0 0.0 0.0 0.0 0.0]

# get the distances all at once, which may or may not be useful...
function readlatlon()
    """ for reasons I don't understand, this can be about 1km different from obspy """
    f=open(LOCATION_FILE)
    readline(f)
    locations=Dict()
    while !eof(f)
        tmp=split(readline(f),",")
        locations[tmp[2]]=(parse(Float32,tmp[3]),parse(Float32,tmp[4]))
    end

    distances=Dict()
    for sta1 in keys(locations)
        for sta2 in keys(locations)
            if sta1==sta2
                distances[sta1*"_"*sta2]=0.
            else
                distances[sta1*"_"*sta2]=Geodesics.surface_distance(
                    locations[sta1][2],locations[sta1][1],
                    locations[sta2][2],locations[sta2][1],6371.
                )
                distances[sta2*"_"*sta1]=distances[sta1*"_"*sta2]
            end
        end
    end
    locations,distances
end
const locations,distances=readlatlon()

function loadHDF5_seg!(fft1,dat,name,work1,work2,iseg,filespaceid,memspaceid)
    #dim0 in those select_hyperslab goes to the right (instead of going downwards) here...
    local nt=nt_wrong

    obj=dat[name*".real"]
    HDF5.h5s_select_hyperslab(filespaceid,HDF5.H5S_SELECT_SET,[iseg-1,0],[1,1],[1,nt],C_NULL)
    HDF5.h5d_read(obj.id,HDF5.hdf5_type_id(Float32),memspaceid,filespaceid,obj.xfer,work1)
    HDF5.refresh(obj)

    obj=dat[name*".imag"]
    HDF5.h5s_select_hyperslab(filespaceid,HDF5.H5S_SELECT_SET,[iseg-1,0],[1,1],[1,nt],C_NULL)
    HDF5.h5d_read(obj.id,HDF5.hdf5_type_id(Float32),memspaceid,filespaceid,obj.xfer,work2)
    HDF5.refresh(obj)

    @. fft1=work1+1im*work2
    return
end

function mov_avg!(A,n,B)
    nlen::Int32=0
    sum::Float32=0.
    n2::Int32=n*2+1
    if (length(A)<n2)
        println("error in mov_avg!!!")
        exit()
    end
    for i=1:n2-1
        sum+=abs(A[i])
    end
    i1::Int32=n2-1
    i2::Int32=1
    for i::Int32=n+1:length(A)-n
        i1+=1
        sum=sum+abs(A[i1])
        B[i]=sum/n2
        sum-=abs(A[i2])
        i2+=1
    end
    for i=1:n
        B[i]=abs(A[i])
    end
    for i=length(A)-n+1:length(A)
        B[i]=abs(A[i])
    end

    #for the stupid last element
    B[end]=1.0

end

function prepare(fft1,fft2,fft1abs,fftcommon,fftplan,ifftplan,n1,n2,c1fft,work1,work2,c1time)
    """ 
    prepare for the c3, which includes
    1. compute cross spectrum, conj(fft1).*conj(fft2)./fft1abs
    2. ifft to time domain
    3. cut the positive part from n1 for nlen samples
    4. fft back to the freq domain
    """
    t1=time()
#    @. work1=conj(fft1)*fft2/fft1abs.^2
#    @. work1=conj(fft1)*fft2/fft1abs/fftcommon
    # @. work1=conj(fft1)*fft2/fft1abs.^2
    @. work1=(conj(fft1)/abs(fft1))*(fft2/abs(fft2))
    if any(isnan,work1) 
        println("stop!")
        exit()
    end
    work1[1]=0.0
    t2=time()
    timeall[1]+=t2-t1

    t1=time()
    FFTW.mul!(work2,ifftplan,work1)
    t2=time()
    timeall[2]+=t2-t1

#    FFTW.mul!(c1fft,fftplan,work2[n1:n2])   #I don't know why this is not working

    t1=time()
    c1fft.=fftplan*work2[n1:n2]
    t2=time()
    timeall[3]+=t2-t1

    c1time.=work2
end


function c3(colpair,colindx,ncol,c1fft,c3stack,n3stack)
    t1=time()
    for i=1:ncol
        tmp1=split(colpair[i],"_")

        for j=i+1:ncol
            tmp2=split(colpair[j],"_")
#            println("tmp2: ",tmp2)
            if (tmp2[1]*tmp2[2] != tmp1[1]*tmp1[2] ) break end   #the way we store C1 allows for early break
            name=join(["source",tmp1[2],tmp1[3],tmp1[4],tmp2[3],tmp2[4]],"_")
            # println("doing:",(tmp1,tmp2,name))
            idx=colindx[name]
            for k=1:n3
                c3stack[k,idx]+=conj(c1fft[k,i])*c1fft[k,j]
            end
            if (any(isnan,c3stack[:,idx]))
                println("NAN!")
                exit()
            end
            n3stack[idx]+=1
        end
    end
    t2=time()
    timeall[4]+=t2-t1
end


function testc3()
    """ the main function that governs how you load the data """
#    allfiles=filter(x->x[end-2:end]==".h5",readdir(FFT_DIR))


    allfiles=filter(x-> x[end-2:end]==".h5" && any(sta->occursin(sta,x),stalist),readdir(FFT_DIR))

    nfiles=length(allfiles)
    println("we will be doing",nfiles,"files")


    # these are for loading the data
    rawfft=ones(ComplexF32,nt,nfiles*length(comp_a))
    work1=Array{Float32}(undef,nt_wrong)
    work2=Array{Float32}(undef,nt_wrong)
    filespaceid=HDF5.h5s_create_simple(2,[nseg,nt_wrong],C_NULL)
    memspaceid=HDF5.h5s_create_simple(1,[nt_wrong],C_NULL)


    # these are for prepare subroutine
    mulfft=Array{ComplexF32}(undef,nt,nfiles^2*length(comp_r)*length(comp_s))
    c1fft=Array{ComplexF32}(undef,div(nlen,2)+1,nfiles^2*length(comp_r)*length(comp_s))    
    absfft=Array{Float32}(undef,nt)
    absfft_common=similar(absfft)
    c1time=Array{Float32}(undef,nt_original,nfiles^2*length(comp_r)*length(comp_s))
    work3=Array{ComplexF32}(undef,nt)
    work4=Array{Float32}(undef,nt_original)
    ifftplan=plan_irfft(work3,nt_original)
    fftplan=plan_rfft(work4[n1:n2])


    # stores the index for the final c3 stack
    colindx_c1=Dict{String,Int32}()
    colindx_c3=Dict{String,Int32}()
    newindx_c1=0
    newindx_c3=0
    colnsta=Array{String}(undef,nfiles*length(comp_a))
    colpair=Array{String}(undef,nfiles^2*length(comp_s)*length(comp_r))
    c1stack=similar(c1time)
    n1stack=zeros(Int32,nfiles^2*length(comp_r)*length(comp_r))
    c3stack=similar(c1fft)
    n3stack=zeros(Int32,nfiles^2*length(comp_s)*length(comp_r)*length(comp_r))
    ifftplan_c3=plan_irfft(c3stack[:,1],nlen)

    h5files=Array{HDF5File}(undef,nfiles)
    namerefs=Array{String}(undef,nfiles)
    datasetnames=Array{Array{String,1},1}()
    for (idx,file) in enumerate(allfiles)
        h5files[idx]=h5open(FFT_DIR*file,"r")
        push!(datasetnames,names(h5files[idx]))
        namerefs[idx]=datasetnames[end][4][1:end-5]
    end


    for day=map(x->Dates.format(x,"yyyy_mm_dd"),d1:Dates.Day(1):d2)
        println("doing day:",day)
        for iseg=1:nseg
            println("doing seg:",iseg)
            # this bits load the raw fft into matrix
            t1=time()
            ndx=0
            for (datasetname,nameref,h5file) in zip(datasetnames,namerefs,h5files)
                tmp=split(nameref,"_")
                for comp in comp_a
                    name=join([tmp[1],tmp[2],tmp[3],comp,day],"_")
                    if name*".real" in datasetname
                        if read(attrs(h5file[name*".real"]),"std")[iseg]>=10
#                            println("std large:",(day,iseg,name))
                            continue
                        end
                        ndx+=1
                        # println("doing:",day," ",name," ",iseg)
                        loadHDF5_seg!(view(rawfft,1:nt_wrong,ndx),h5file,name,work1,work2,iseg,filespaceid,memspaceid)
                        if any(iszero,view(rawfft,1:nt_wrong,ndx))
                            println(("zeros!!!!",day,name,iseg))
                            ndx-=1
                            continue
                        end
                        colnsta[ndx]=tmp[3]*"_"*comp
                    end
                end
            end
            t2=time()
            timeall[5]+=t2-t1

#            if ndx<=3 continue end
            if ndx<=3 continue end

            # now prepare the c1; also establish the index if this pair doesn't exist before
            ncol=0
            iassigned=false
            for i=1:ndx
                if !(split(colnsta[i],"_")[2] in comp_s) continue end
                mov_avg!(view(rawfft,:,i),10,absfft)
                if (iassigned==false)
                    absfft_common.=absfft
                    iassigned=true
                end
                for j=1:ndx
                    if !(split(colnsta[j],"_")[2] in comp_r) continue end
                    ncol+=1

                    prepare(view(rawfft,:,i),view(rawfft,:,j),absfft,absfft_common,fftplan,ifftplan,
                        n1,n2,view(c1fft,:,ncol),work3,work4,view(c1time,:,ncol))

                    col=colnsta[i]*"_"*colnsta[j]
                    colpair[ncol]=col
                end
            end


            # set up the index for C1 and C3 independently
            for i=1:ndx
                if !(split(colnsta[i],"_")[2] in comp_r) continue end
                for j=1:ndx
                    if !(split(colnsta[j],"_")[2] in comp_r) continue end

                    col=colnsta[i]*"_"*colnsta[j]
                    if (col in keys(colindx_c1)) == false
                        newindx_c1+=1
                        colindx_c1[col]=newindx_c1
                        println(("adding c1:",newindx_c1,col))
                    end

                    for comp in comp_s
                        col="source_"*comp*"_"*colnsta[i]*"_"*colnsta[j]
                        if (col in keys(colindx_c3)) == false
                            newindx_c3+=1
                            colindx_c3[col]=newindx_c3
                            println(("adding c3:",newindx_c3,col))
                        end
                    end
                end
            end

            # let's do c3
            c3(colpair,colindx_c3,ncol,c1fft,c3stack,n3stack)

        end   # end the loop over segment

        for key in keys(colindx_c3)
            idx=colindx_c3[key]
            tmp=split(key,"_")
            tr=SAC.sample()
            tr.delta=1.0/downsamp_freq
            tr.evla=locations[tmp[3]][1]
            tr.evlo=locations[tmp[3]][2]
            tr.stla=locations[tmp[5]][1]
            tr.stlo=locations[tmp[5]][2]
            tr.kstnm=tmp[3]*tmp[5]
            if (n3stack[idx]>0)
                tr.npts=nlen
                println(any(isnan,c3stack[:,idx]))
                c3stack[:,idx]/=n3stack[idx]
                tr.t=fftshift(ifftplan_c3*c3stack[:,idx])
                write(tr,CCF_DIR*key*"_"*day*".C3.SAC",byteswap=false)
            end
        end
        c3stack.=0.
        n3stack.=0       
        c1fft.=0.
        n1stack.=0

    end
            
    # close all h5files and ends
    foreach(close,h5files)
    HDF5.h5s_close(memspaceid)
    HDF5.h5s_close(filespaceid)
end

testc3()

println(timeall)
