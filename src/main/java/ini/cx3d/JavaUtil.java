package ini.cx3d;

import ini.cx3d.spatialOrganization.NewDelaunayTest;
import ini.cx3d.spatialOrganization.OpenTriangleOrganizer;
import ini.cx3d.swig.spatialOrganization.JavaUtilT_PhysicalNode;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * provides functionality that has not been implemented yet in C++
 * especially static methods as they can't be handled by SWIG directors
 */
public class JavaUtil extends JavaUtilT_PhysicalNode {

    @Override
    public OpenTriangleOrganizer oto_createSimpleOpenTriangleOrganizer() {
        return OpenTriangleOrganizer.createSimpleOpenTriangleOrganizer_java();
    }

    static List<Integer> list = Arrays.asList(new Integer[]{
            0, 1, 2, 3});

    @Override
    public int[] generateTriangleOrder() {
        Collections.shuffle(list, NewDelaunayTest.rand);
        int[] ret = new int[4];
        for (int i = 0; i < list.size(); i++) {
            ret[i] = list.get(i);
        }
        return ret;
    }
}
