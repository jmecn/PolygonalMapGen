package net.jmecn;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import net.jmecn.voronoi.Center;
import net.jmecn.voronoi.Edge;
import net.jmecn.voronoi.Voronoi;

public class Main {

    public static void main(String[] args) {
        Voronoi v = new Voronoi();
        v.InitVoroni(500, 513, 513);

        BufferedImage image = new BufferedImage(400, 400, BufferedImage.TYPE_4BYTE_ABGR);

        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, 400, 400);

        g.setColor(new Color(0, 0, 0));
        for (Edge e : v.voronoiEdgeList) {
            g.drawLine((int)e.a.x, (int)e.a.y, (int)e.b.x, (int)e.b.y);
        }
        
        try {
            ImageIO.write(image, "png", new File("1.png"));
        } catch (IOException e1) {
            e1.printStackTrace();
        }

        g.setColor(Color.RED);
        for (Center c : v.voronoiCenterList) {
            g.drawOval((int)c.point.x-1, (int)c.point.y-1, 2, 2);
        }
        
        try {
            ImageIO.write(image, "png", new File("2.png"));
        } catch (IOException e1) {
            e1.printStackTrace();
        }


    }

}
