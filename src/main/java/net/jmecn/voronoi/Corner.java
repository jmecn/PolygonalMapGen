package net.jmecn.voronoi;

import java.util.HashMap;
import java.util.Map;

public class Corner {
    public int index;

    public Site point; // location
    public boolean water; // lake or ocean
    public boolean ocean; // ocean
    public boolean coast; // land polygon touching an ocean
    public boolean border; // at the edge of the map
    public float elevation = 0; // 0.0-1.0
    public float moisture = 0; // 0.0-1.0

    public Map<Integer, Center> touches;
    public Map<Edge, Edge> protrudes;
    public Map<Integer, Corner> adjacent;

    public int river; // 0 if no river, or volume of water in river
    public Corner downslope; // pointer to adjacent corner most downhill
    public Corner watershed; // pointer to coastal corner, or null
    public int watershed_size;
    public int water_distance;

    public Corner() {
        touches = new HashMap<Integer, Center>();
        protrudes = new HashMap<Edge, Edge>();
        adjacent = new HashMap<Integer, Corner>();
    }
}
