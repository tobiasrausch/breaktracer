#ifndef CLUSTER_H
#define CLUSTER_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace breaktracer
{

  // Edge struct
  template<typename TWeight, typename TVertex>
  struct EdgeRecord {
    typedef TVertex TVertexType;
    TVertex source;
    TVertex target;
    TWeight weight;
    
    EdgeRecord(TVertex s, TVertex t, TWeight w) : source(s), target(t), weight(w) {}

    bool operator<(const EdgeRecord e2) const {
      return ((weight < e2.weight) || ((weight == e2.weight) && (source < e2.source)) || ((weight == e2.weight) && (source == e2.source) && (target < e2.target)));
    }
  };


  template<typename TConfig, typename TCompEdgeList>
  inline void
  _searchCliques(TConfig const& c, TCompEdgeList& compEdge, std::vector<TraceRecord>& tr, std::vector<BrInTrace>& sv) {
    typedef typename TCompEdgeList::mapped_type TEdgeList;
    typedef typename TEdgeList::value_type TEdgeRecord;
    typedef typename TEdgeRecord::TVertexType TVertex;

    // Iterate all components
    for(typename TCompEdgeList::iterator compIt = compEdge.begin(); compIt != compEdge.end(); ++compIt) {
      // Sort edges by weight
      std::sort(compIt->second.begin(), compIt->second.end());

      // Find a large clique
      typename TEdgeList::const_iterator itWEdge = compIt->second.begin();
      typename TEdgeList::const_iterator itWEdgeEnd = compIt->second.end();
      typedef std::set<TVertex> TCliqueMembers;
      typedef std::set<std::size_t> TSeeds;
      TCliqueMembers clique;
      TCliqueMembers incompatible;
      TSeeds seeds;
      
      // Initialize clique
      clique.insert(itWEdge->source);
      seeds.insert(tr[itWEdge->source].id);
      int32_t chr = tr[itWEdge->source].chr;
      int32_t chr2 = tr[itWEdge->source].chr2;
      int32_t ciposlow = tr[itWEdge->source].pos;
      std::vector<int32_t> refpos;
      refpos.push_back(tr[itWEdge->source].pos);
      int32_t ciposhigh = tr[itWEdge->source].pos; 
      int32_t ciendlow = tr[itWEdge->source].pos2;
      std::vector<int32_t> refpos2;
      refpos2.push_back(tr[itWEdge->source].pos2);
      int32_t ciendhigh = tr[itWEdge->source].pos2;
      std::vector<int32_t> qual;
      qual.push_back(tr[itWEdge->source].qual);
      std::vector<int32_t> inslen;
      inslen.push_back(tr[itWEdge->source].inslen);

      // Initialize wiggle
      uint32_t wiggle = c.minRefSep;

      // Grow clique
      bool cliqueGrow = true;
      while (cliqueGrow) {
	itWEdge = compIt->second.begin();
	cliqueGrow = false;
	// Find next best edge for extension
	for(;(!cliqueGrow) && (itWEdge != itWEdgeEnd);++itWEdge) {
	  TVertex v;
	  if ((clique.find(itWEdge->source) == clique.end()) && (clique.find(itWEdge->target) != clique.end())) v = itWEdge->source;
	  else if ((clique.find(itWEdge->source) != clique.end()) && (clique.find(itWEdge->target) == clique.end())) v = itWEdge->target;
	  else continue;
	  if (incompatible.find(v) != incompatible.end()) continue;
	  if (seeds.find(tr[v].id) != seeds.end()) continue;
	  // Try to update clique with this vertex
	  int32_t newCiPosLow = std::min(tr[v].pos, ciposlow);
	  int32_t newCiPosHigh = std::max(tr[v].pos, ciposhigh);
	  int32_t newCiEndLow = std::min(tr[v].pos2, ciendlow);
	  int32_t newCiEndHigh = std::max(tr[v].pos2, ciendhigh);
	  if (((newCiPosHigh - newCiPosLow) < (int32_t) wiggle) && ((newCiEndHigh - newCiEndLow) < (int32_t) wiggle)) cliqueGrow = true;
	  if (cliqueGrow) {
	    // Accept new vertex
	    clique.insert(v);
	    if (seeds.find(tr[v].id) == seeds.end()) {
	      seeds.insert(tr[v].id);
	      ciposlow = newCiPosLow;
	      refpos.push_back(tr[v].pos);
	      ciposhigh = newCiPosHigh;
	      ciendlow = newCiEndLow;
	      refpos2.push_back(tr[v].pos2);
	      ciendhigh = newCiEndHigh;
	      qual.push_back(tr[v].qual);
	      inslen.push_back(tr[v].inslen);
	    }
	  } else incompatible.insert(v);
	}
      }

      // New SV-INS-SV?
      if (seeds.size() >= c.minCliqueSize) {
	int32_t fpos = mdv(refpos);
	int32_t fpos2 = mdv(refpos2);
	int32_t qout = mdv(qual);
	int32_t iout = mdv(inslen);
	sv.push_back(BrInTrace(chr, fpos, chr2, fpos2, qout, iout, seeds));
      }
    }
  }
  

  template<typename TConfig>
  inline void
  cluster(TConfig const& c, std::vector<TraceRecord>& tr, std::vector<BrInTrace>& sv) {
    uint32_t count = 0;
    for(int32_t refIdx = 0; refIdx < c.nchr; ++refIdx) {
      
      // Components
      typedef std::vector<uint32_t> TComponent;
      TComponent comp;
      comp.resize(tr.size(), 0);
      uint32_t numComp = 0;

      // Edge lists for each component
      typedef uint32_t TWeightType;
      typedef uint32_t TVertex;
      typedef EdgeRecord<TWeightType, TVertex> TEdgeRecord;
      typedef std::vector<TEdgeRecord> TEdgeList;
      typedef std::map<uint32_t, TEdgeList> TCompEdgeList;
      TCompEdgeList compEdge;

	
      std::size_t lastConnectedNode = 0;
      std::size_t lastConnectedNodeStart = 0;
      for(uint32_t i = 0; i<tr.size(); ++i) {
	if (tr[i].chr == refIdx) {
	  ++count;
	  // Safe to clean the graph?
	  if (i > lastConnectedNode) {
	    // Clean edge lists
	    if (!compEdge.empty()) {
	      // Search cliques
	      _searchCliques(c, compEdge, tr, sv);
	      lastConnectedNodeStart = lastConnectedNode;
	      compEdge.clear();
	    }
	  }

	  // Get size variability
	  for(uint32_t j = i + 1; j<tr.size(); ++j) {
	    if (tr[j].chr == refIdx) {
	      if ( std::abs(tr[j].pos - tr[i].pos) > c.minRefSep) break;
	      if ( std::abs(tr[j].pos2 - tr[i].pos2) > c.minRefSep) break;
	      if ( std::abs(tr[j].inslen - tr[i].inslen) <= c.minRefSep) {
		// Update last connected node
		if (j > lastConnectedNode) lastConnectedNode = j;
		
		// Assign components
		uint32_t compIndex = 0;
		if (!comp[i]) {
		  if (!comp[j]) {
		    // Both vertices have no component
		    compIndex = ++numComp;
		    comp[i] = compIndex;
		    comp[j] = compIndex;
		    compEdge.insert(std::make_pair(compIndex, TEdgeList()));
		  } else {
		    compIndex = comp[j];
		    comp[i] = compIndex;
		  }	
		} else {
		  if (!comp[j]) {
		    compIndex = comp[i];
		    comp[j] = compIndex;
		  } else {
		    // Both vertices have a component
		    if (comp[j] == comp[i]) {
		      compIndex = comp[j];
		    } else {
		      // Merge components
		      compIndex = comp[i];
		      uint32_t otherIndex = comp[j];
		      if (otherIndex < compIndex) {
			compIndex = comp[j];
			otherIndex = comp[i];
		      }
		      // Re-label other index
		      for(uint32_t k = lastConnectedNodeStart; k <= lastConnectedNode; ++k) {
			if (otherIndex == comp[k]) comp[k] = compIndex;
		      }
		      // Merge edge lists
		      TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
		      TCompEdgeList::iterator compEdgeOtherIt = compEdge.find(otherIndex);
		      compEdgeIt->second.insert(compEdgeIt->second.end(), compEdgeOtherIt->second.begin(), compEdgeOtherIt->second.end());
		      compEdge.erase(compEdgeOtherIt);
		    }
		  }
		}
		
		// Append new edge
		TCompEdgeList::iterator compEdgeIt = compEdge.find(compIndex);
		if (compEdgeIt->second.size() < c.graphPruning) {
		  // SV-INS-SV distance
		  TWeightType weight = std::abs(tr[j].pos2 - tr[i].pos2) + std::abs(tr[j].pos - tr[i].pos) + std::abs(tr[j].inslen - tr[i].inslen);
		  compEdgeIt->second.push_back(TEdgeRecord(i, j, weight));
		}
	      }
	    }
	  }
	}
      }
      // Search cliques
      if (!compEdge.empty()) {
	_searchCliques(c, compEdge, tr, sv);
	compEdge.clear();
      }
    }
  }


}

#endif
