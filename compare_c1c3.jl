using Dates
using SAC
using Glob
using Statistics
using PlotlyJS
using ORCA

d1=Date(2012,1,1)
d2=Date(2012,12,30)
yr="2012"
CURRDIR="/Users/zma/v2.0_4Xperformace"

# remember the prefix may have _ at the end and postfix may have . at the beginning
function cal_corr(xmin,xmax,pair,WORKDIR,prefix_ref,postfix_ref,prefix,postfix,f1,f2)
    cd(WORKDIR)
    file=SAC.read(prefix_ref*pair*postfix_ref)
    x=range(-file.npts*file.delta/2.0,stop=file.npts*file.delta/2.0,length=file.npts)
    idx=(x.>xmin) .& (x.<xmax)
    yref=SAC.bp(file,f1,f2).t[idx]
    xplot=Array{String,1}()
    yplot=Array{Float32,1}()

    for day=map(x->Dates.format(x,"yyyy_mm_dd"),d1:Dates.Day(1):d2)
       file=prefix*pair*"_"*day*postfix
       if !(isfile(file)) continue end
       y=SAC.bp(SAC.read(file),f1,f2).t[idx]
       push!(xplot,day[6:10])
       push!(yplot,simple_corr(y,yref))
    end
    cd(CURRDIR)
    xplot,yplot
end

function simple_corr(x,y)
    return sum(x.*y)/sqrt(sum(x.*x)*sum(y.*y))
end

function cal_sum_corr(xmin,xmax,pair,WORKDIR,prefix_ref,postfix_ref,prefix,postfix,f1,f2)
    cd(WORKDIR)
    file=SAC.read(prefix_ref*pair*postfix_ref)
    x=range(-file.npts*file.delta/2.0,stop=file.npts*file.delta/2.0,length=file.npts)
    idx=(x.>xmin) .& (x.<xmax)
    yref=SAC.bp(file,f1,f2).t[idx]
    xplot=Array{String,1}()
    yplot=Array{Float32,1}()
    trace_sum=zeros(Float32,size(yref))

    for day=map(x->Dates.format(x,"yyyy_mm_dd"),d1:Dates.Day(1):d2)
       file=prefix*pair*"_"*day*postfix
       if !(isfile(file)) continue end
       y=SAC.bp(SAC.read(file),f1,f2).t[idx]
       @. trace_sum+=y
       push!(xplot,day[6:10])
       push!(yplot,simple_corr(trace_sum,yref))
    end
    cd(CURRDIR)
    xplot,yplot
end


# simply extract the actual cross correlation function
function extract_corr(xmin,xmax,pair,WORKDIR,prefix_ref,postfix_ref,f1,f2)
    cd(WORKDIR)
    file=SAC.read(prefix_ref*pair*postfix_ref)
    x=range(-file.npts*file.delta/2.0,stop=file.npts*file.delta/2.0,length=file.npts)
    idx=(x.>xmin) .& (x.<xmax)
    y=SAC.bp(file,f1,f2).t[idx]
    y/=maximum(y)
    cd(CURRDIR)
    return x[idx],y
end

function getpairs(WORKDIR)
    cd(WORKDIR)
    pairs=Array{String,1}()
    for file in filter(x->x[1:6]!="source",glob("*stack.SAC"))
        pair=split(file,".")[1]
        tmp=split(pair,"_")
        if (tmp[1]>=tmp[3]) continue end
        if !(pair in pairs) push!(pairs,pair) end
    end
    pairs 
end

# make some plots for all pairs
function plot_allpairs()
    xmin=0.
    xmax=100.
    f1=0.1
    f2=0.3

    WORKDIR = "/Users/zma/v2.0_4Xperformace/testing.dir/CCF"
    mypairs=getpairs(WORKDIR)
    rr_c1=zeros((d2-d1).value+1,length(mypairs)) 
    for (iy,pair) in enumerate(mypairs)
        println(("doing:",pair))
        xplots,yplots=cal_sum_corr(xmin,xmax,pair,WORKDIR,"",".C1stack.SAC","",".C1.SAC",f1,f2)
        for (xplot,yplot) in zip(xplots,yplots)
            rr_c1[(Date("2012_"*xplot,"yy_mm_dd")-d1).value+1,iy]=yplot
        end
        fillgap(view(rr_c1,:,iy),0.0)
    end

#    WORKDIR = "/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/"
    WORKDIR = "/Users/zma/v2.0_4Xperformace/testing.dir/n150nlen1200/"
    rr_c3=zeros(size(rr_c1))
    for (iy,pair) in enumerate(mypairs)
        println(("doing:",pair))
        xplots,yplots=cal_sum_corr(xmin,xmax,pair,WORKDIR,"source_BHZ_",".C3stack.SAC",
            "source_BHZ_",".C3.SAC",f1,f2)
        for (xplot,yplot) in zip(xplots,yplots)
            rr_c3[(Date("2012_"*xplot,"yy_mm_dd")-d1).value+1,iy]=yplot
        end
        fillgap(view(rr_c3,:,iy),0.0)
    end

    p1 = plot(heatmap(;y=mypairs,z=rr_c1),Layout(;width=500,height=800,margin_l=150,title="c1"))
    p2 = plot(heatmap(;y=mypairs,z=rr_c3),Layout(;width=500,height=800,margin_l=150,title="c3"))
    p3 = plot(heatmap(;y=mypairs,z=rr_c3-rr_c1),Layout(;width=500,height=800,margin_l=150,title="c3-c1"))
    savefig(p1,"allpairs_c1_0.1-0.3_n150nlen1200.pdf")
    savefig(p2,"allpairs_c3_0.1-0.3_n150nlen1200.pdf")
    savefig(p3,"allpairs_c3-c1_0.1-0.3_n150nlen1200.pdf")
    pp=display([p1,p2,p3])
    #return mypairs,rr_c1,rr_c3
end

# for plotting heatmap, fill the gap in vector a
# assume there's at least one nongap element
# i point to 0 or the last nongap element
# j point to gap
function fillgap(a,gap)
    j0=findfirst(isequal(gap),a)
    if (j0==nothing) return end
    i=j0-1
    for j=j0:length(a)+1
        if (j>length(a)) || (a[j]!=gap)
           if (j>length(a))
               a[i+1:end].=a[i]
           else  
               if (i==0)
                   a[i+1:j-1].=a[j]
               else
                   k=div(i+j,2)
                   a[i+1:k].=a[i]
                   a[k+1:j-1].=a[j]
               end
               i=j
           end
        end
    end
end


# this part actually makes the plots for publication
# the daily plot uses plot_sac_daily.jl
function main()
    pair="GNW_BHZ_STOR_BHZ"
    f1=0.1
    f2=0.3

    # plot the actual cross correlation
    xmin=-200.
    xmax=200.
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/CCF"
    xplot_c1,yplot_c1=extract_corr(xmin,xmax,pair,WORKDIR,"",".C1stack.SAC",f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/"
    xplot_c3_10_BHZ,yplot_c3_10_BHZ=extract_corr(xmin,xmax,pair,WORKDIR,"source_BHZ_",".C3stack.SAC",
         f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/"
    xplot_c3_10,yplot_c3_10=extract_corr(xmin,xmax,pair,WORKDIR,"",".C3stack.SAC",
         f1,f2)
    tr1=scatter(;x=xplot_c1,y=yplot_c1,mode="lines",name="C1")
    tr2=scatter(;x=xplot_c3_10_BHZ,y=yplot_c3_10_BHZ,mode="lines",name="C3 10 sources BHZ only")   
    tr3=scatter(;x=xplot_c3_10,y=yplot_c3_10,mode="lines",name="C3 10 sources")   
    pp=plot([tr2;tr3;tr1],Layout(;width=700,height=300,autosize=false))
    savefig(pp,"CC.pdf")

    idx = (xplot_c1.>0) .& (xplot_c1.<100)
    println(simple_corr(yplot_c1[idx],yplot_c3_10_BHZ[idx]))
    println(simple_corr(yplot_c1[idx],yplot_c3_10[idx]))


    # plot the cal_corr
    xmin=0.
    xmax=100.
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/CCF"
    xplot_c1,yplot_c1=cal_corr(xmin,xmax,pair,WORKDIR,"",".C1stack.SAC","",".C1.SAC",f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/" 
    xplot_c3_10_BHZ,yplot_c3_10_BHZ=cal_corr(xmin,xmax,pair,WORKDIR,"source_BHZ_",".C3stack.SAC",
          "source_BHZ_",".C3.SAC",f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/" 
    xplot_c3_10,yplot_c3_10=cal_corr(xmin,xmax,pair,WORKDIR,"",".C3stack.SAC",
          "source_BHZ_",".C3.SAC",f1,f2)
    tr1=scatter(;x=xplot_c1,y=yplot_c1,mode="lines",name="C1")
    tr2=scatter(;x=xplot_c3_10_BHZ,y=yplot_c3_10_BHZ,mode="lines",name="C3 10 sources BHZ only")   
    tr3=scatter(;x=xplot_c3_10,y=yplot_c3_10,mode="lines",name="C3 10 sources")   
    pp=plot([tr2;tr3;tr1],Layout(;width=600,height=700,autosize=false))

    # plot the cal_sum_corr
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/CCF"
    xplot_c1,yplot_c1=cal_sum_corr(xmin,xmax,pair,WORKDIR,"",".C1stack.SAC","",".C1.SAC",f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/"
    xplot_c3_10_BHZ,yplot_c3_10_BHZ=cal_sum_corr(xmin,xmax,pair,WORKDIR,"source_BHZ_",".C3stack.SAC",
          "source_BHZ_",".C3.SAC",f1,f2)
    WORKDIR="/Users/zma/v2.0_4Xperformace/testing.dir/10source_BHENZ_daily/"
    xplot_c3_10,yplot_c3_10=cal_sum_corr(xmin,xmax,pair,WORKDIR,"",".C3stack.SAC",
          "source_BHZ_",".C3.SAC",f1,f2)
    tr1=scatter(;x=xplot_c1,y=yplot_c1,mode="lines",name="C1")
    tr2=scatter(;x=xplot_c3_10_BHZ,y=yplot_c3_10_BHZ,mode="lines",name="C3 10 sources BHZ only")   
    tr3=scatter(;x=xplot_c3_10,y=yplot_c3_10,mode="lines",name="C3 10 sources")   
    pp=plot([tr2;tr3;tr1],Layout(;width=600,height=700,autosize=false))
    savefig(pp,"sumcorr.pdf")
end

main()
