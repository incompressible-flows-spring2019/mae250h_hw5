export regularize!, interpolate!

roma1(r) = (1+sqrt(-3r^2+1))/3
roma2(r) = (5-3r-sqrt(1-3(1-r)^2))/6

@inline ddf_roma(r::Real) = r > 1.5 ? 0.0 : r <= 0.5 ? roma1(r) : roma2(r)

@inline ddf_test(r::Real) = r^2


"""
    regularize!(s::NodeData,f::ScalarPoint,pts::VectorPoint,p::Parameters)

Regularize the scalar data `f` at points `pts` to nodes on the grid whose parameters
are provided in `p`. Return the resulting regularized field in `s`.
"""
function regularize!(s::NodeData{NX,NY},f::ScalarPoint{N},pts::VectorPoint{N},p::Parameters) where {NX,NY,N}

  i_n = indices(s,1,interior=false)
  j_n = indices(s,2,interior=false)

  xg = xmap(i_n,s,p)
  yg = ymap(j_n,s,p)
  scale = 1/p.Δx

  for (j,yj) in enumerate(yg)
    for (i,xi) in enumerate(xg)
      stmp = 0.0
      for pt in 1:N
        stmp += f[pt]*ddf_roma(abs(pts.x[pt]-xi)*scale)*ddf_roma(abs(pts.y[pt]-yj)*scale)
      end
      s[i,j] += stmp
    end
  end
  return s
end

"""
    interpolate!(f::ScalarPoint,s::NodeData,pts::VectorPoint,p::Parameters)

Interpolate the scalar node data `s` on the grid whose parameters
are provided in `p`, to points `pts`. Return the resulting
interpolated data in `f`.
"""
function interpolate!(f::ScalarPoint{N},s::NodeData{NX,NY},pts::VectorPoint{N},p::Parameters) where {NX,NY,N}

  i_n = indices(s,1,interior=false)
  j_n = indices(s,2,interior=false)

  xg = xmap(i_n,s,p)
  yg = ymap(j_n,s,p)
  scale = 1/p.Δx

  for (j,yj) in enumerate(yg)
    for (i,xi) in enumerate(xg)
      stmp = s[i,j]
      for pt in 1:N
        f[pt] += stmp*ddf_roma(abs(pts.x[pt]-xi)*scale)*ddf_roma(abs(pts.y[pt]-yj)*scale)
      end
    end
  end
  return f
end
