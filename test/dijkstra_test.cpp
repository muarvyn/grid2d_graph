#include <iostream>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#define NEED_PROPERTY_MAP
#include "grid2d_graph.hpp"

const int ROWS_NUM=4, COLS_NUM=5;
typedef g2dg::grid2d_graph<int, int, ROWS_NUM, COLS_NUM> graph_t;

int main (int argc, char const *argv[]) {
  using namespace boost;
  // Check the concepts that graph models.  This is included to demonstrate
  // how concept checking works, but is not required for a working program
  // since Boost algorithms do their own concept checking.
  #ifdef NEED_IN_EDGE_ITERATOR
  BOOST_CONCEPT_ASSERT(( BidirectionalGraphConcept<graph_t> ));
  #endif
  #ifdef NEED_ADJACENCY_ITERATOR
  BOOST_CONCEPT_ASSERT(( AdjacencyGraphConcept<graph_t> ));
  #endif
  BOOST_CONCEPT_ASSERT(( VertexListGraphConcept<graph_t> ));
#ifdef NEED_EDGE_ITERATOR
  BOOST_CONCEPT_ASSERT(( EdgeListGraphConcept<graph_t> ));
#endif
#ifdef NEED_ADJACENCY_MATRIX
  BOOST_CONCEPT_ASSERT(( AdjacencyMatrixConcept<graph_t> ));
#endif
#ifdef NEED_PROPERTY_MAP
  BOOST_CONCEPT_ASSERT((
    ReadablePropertyMapConcept<g2dg::edge_weight_map_traits<graph_t>::const_type,
        graph_t::edge_descriptor> ));
#endif
#ifdef NEED_PROPERTY_GRAPH
  BOOST_CONCEPT_ASSERT((
    ReadablePropertyGraphConcept<graph_t, graph_t::edge_descriptor, edge_weight_t> ));
#endif

  graph_t g;

  // Print the outgoing edges of all the vertices.
  //
  std::cout << "Vertices, outgoing edges, and adjacent vertices" << std::endl;
  typedef std::pair< g2dg::g2dg_vertex_iterator<int,int,ROWS_NUM,COLS_NUM>,
      g2dg::g2dg_vertex_iterator<int,int,ROWS_NUM,COLS_NUM> >
      vertex_range;
  //g2dg_vertex_iterator<int,int,ROWS_NUM,COLS_NUM> v_begin, v_end;
  vertex_range range = g2dg::vertices<graph_t>(g);
  graph_t::vertex_iterator vi=range.first, vi_end=range.second;
//  for (boost::tie(vi, vi_end) = g2dg::vertices<ROWS_NUM,COLS_NUM>(g); vi != vi_end; vi++) {
  for (; vi != vi_end; vi++) {
      graph_t::vertex_descriptor u = *vi;
      std::cout << "Vertex " << u << ": ";

      // Adjacenct edges
      typedef std::pair< g2dg::g2dg_out_edge_iterator<int,int,ROWS_NUM,COLS_NUM>,
          g2dg::g2dg_out_edge_iterator<int,int,ROWS_NUM,COLS_NUM> >
          out_edge_range;
      out_edge_range e_range = out_edges(u, g);
      graph_t::out_edge_iterator ei=e_range.first, ei_end=e_range.second;
      // for (boost::tie(ei, ei_end) = g2dg::out_edges<ROWS_NUM,COLS_NUM>(u, g); ei != ei_end; ei++)
      for (; ei != ei_end; ei++)
          std::cout << *ei << "  ";
#ifdef NEED_ADJACENCY_ITERATOR
    std::cout << " Adjacent vertices ";
    // Adjacent vertices
    // Here we want our adjacency_iterator and not boost::adjacency_iterator.
    graph_t::adjacency_iterator ai, ai_end;
    for (boost::tie(ai, ai_end) = ring_graph_ns::adjacent_vertices(u, g); ai != ai_end; ai++) {
      std::cout << *ai << " ";
    }
#endif
    std::cout << std::endl;
  }
  std::cout << g2dg::num_vertices(g) << " vertices" << std::endl << std::endl;

  // Print all the edges in the graph along with their weights.
  //
  #ifdef NEED_EDGE_ITERATOR
  std::cout << "Edges and weights" << std::endl;
  graph_t::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = g2dg::edges(g); ei != ei_end; ei++) {
    graph_t::edge_descriptor e = *ei;
    std::cout << e << " weight " << get(edge_weight, g, e) << std::endl;
  }
  std::cout << g2dg::num_edges(g) << " edges"  << std::endl;
  #endif

    std::cout << std::endl;

    // Do a Dijkstra search from vertex 0.  For n=5 this will print:
    //

    graph_t::vertex_descriptor s(0,0);
    /*
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    for (std::size_t j = 0; j < ?; ++j) {
        weightmap[e] = weights[j];
    }*/

    std::vector<int> wmap_init=
        {1, 1, 0, 0, 1,
         1, 1, 1, 0, 1,
         1, 1, 1, 0, 0,
         0, 0, 0, 1, 0
        };
    g2dg::edge_weight_map<graph_t> wmap(wmap_init.begin(),wmap_init.end());
    std::vector<graph_t::vertex_descriptor> pred(num_vertices(g));
    std::vector<g2dg::edge_weight_map_traits<graph_t>::value_type> dist(num_vertices(g));
    iterator_property_map<std::vector<graph_t::vertex_descriptor>::iterator,
                          property_map<graph_t, vertex_index_t>::const_type>
        pred_pm(pred.begin(), get(vertex_index, g));
    iterator_property_map<std::vector<
        g2dg::edge_weight_map_traits<graph_t>::value_type>::iterator,
        property_map<graph_t, vertex_index_t>::const_type
    >
        dist_pm(dist.begin(), get(vertex_index, g));

    dijkstra_shortest_paths(g, s,
                            predecessor_map(pred_pm).
                            distance_map(dist_pm).
                            weight_map(wmap) );

    property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);
    std::cout << "Dijkstra search from vertex " << s << std::endl;
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      graph_t::vertex_descriptor u = *vi;
      std::cout << "Vertex " << u << ": "
                << "parent "<< pred[indexmap[u]] << ", "
                << "distance " << dist[indexmap[u]]
                << std::endl;
    }

    return 0;
}
