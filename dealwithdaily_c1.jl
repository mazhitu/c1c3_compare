using Dates
using SAC
using Glob
using Statistics

d1=Date(2012,1,1)
d2=Date(2012,12,30)

WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/CCF/"
CURRDIR="/Users/zma/v2.0_4Xperformace"

function remove_large_std()

    cd(WORKDIR)

    pairs=Array{String,1}()
    for file in glob("*.C1.SAC")
        tmp=split(file,"_")
        pair=join(tmp[1:4],"_")
        if !(pair in pairs) push!(pairs,pair) end
    end
    println("number of pairs:",length(pairs))

    for pair in pairs
        stds=Array{Float32,1}()
        days=Array{String,1}()
        for day=map(x->Dates.format(x,"yyyy_mm_dd"),d1:Dates.Day(1):d2)
            file=pair*"_"*day*".C1.SAC"
            if !(isfile(file)) continue end
            trace=SAC.read(file)
            push!(stds,std(trace.t))
            push!(days,file)
        end
        if length(stds)<=1 continue end
        # display(plot(days,stds))
        cutoff=median(stds)*10
        println(cutoff)
        
        for (file,std) in zip(days,stds)
            if std>cutoff
                fileout=replace(file,"C1."=>"C1bad.")
                println(file,fileout)
                mv(file,fileout,force=true)
            end
        end
    end

    cd(CURRDIR)
end

function stack_daily()
    cd(WORKDIR)

    pairs=Array{String,1}()
    # for file in filter(x->x[end-5:end]=="C1.SAC",readdir(WORKDIR))
    for file in glob("*.C1.SAC")
        tmp=split(file,"_")
        pair=join(tmp[1:4],"_")
        if !(pair in pairs) push!(pairs,pair) end
    end
    println("number of pairs:",length(pairs))

    for pair in pairs
        trace_sum=nothing
        n_sum=0
        for day=map(x->Dates.format(x,"yyyy_mm_dd"),d1:Dates.Day(1):d2)
            file=pair*"_"*day*".C1.SAC"
            if !(isfile(file)) continue end
            trace=SAC.read(file)
            if (trace_sum==nothing)
                trace_sum=SAC.sample()
                trace_sum.delta=trace.delta
                trace_sum.evla=trace.evla
                trace_sum.evlo=trace.evlo
                trace_sum.stla=trace.stla
                trace_sum.stlo=trace.stlo
                trace_sum.kstnm=trace.kstnm
                trace_sum.npts=trace.npts
                trace_sum.t=trace.t[:]
                n_sum+=1
                println(typeof(trace_sum))
            else
                trace_sum.t=trace_sum.t+trace.t
                n_sum+=1
            end
        end
        trace_sum.t/=n_sum
        tmpfile=pair*".C1stack.SAC"
        if (isfile(tmpfile)) rm(tmpfile,force=true) end
        SAC.write(trace_sum,pair*".C1stack.SAC")
    end

    cd(CURRDIR)
end

remove_large_std()


