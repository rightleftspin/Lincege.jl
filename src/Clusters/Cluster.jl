struct Cluster <: AbstractCluster
    evs::ExpansionVertices
    lc::Rational
    ghash::UInt64
end

function Cluster(evs::ExpansionVertices, cs::TranslationClusters)

end

function Cluster(evs::ExpansionVertices, cs::IsomorphicClusters)

end

Base.length(c::Cluster) = length(c.evs)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
