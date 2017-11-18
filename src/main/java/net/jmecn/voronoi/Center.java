package net.jmecn.voronoi;

import java.util.HashMap;
import java.util.Map;

public class Center {
    public int index;

    public Site point; // location
    public boolean water; // lake or ocean
    public boolean ocean; // ocean
    public boolean coast; // land polygon touching an ocean
    public boolean border; // at the edge of the map
    public BiomeType biome; // biome type (see article)
    public float elevation = 0; // 0.0-1.0
    public float moisture = 0; // 0.0-1.0

    public Map<Integer, Center> neighbors;
    public Map<Edge, Edge> borders;
    public Map<Site, Corner> corners;

    public Center() {
        neighbors = new HashMap<Integer, Center>();
        borders = new HashMap<Edge, Edge>();
        corners = new HashMap<Site, Corner>();
    }
}
