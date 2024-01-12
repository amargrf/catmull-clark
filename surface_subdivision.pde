import controlP5.*;
import peasy.*;

int numIterations = 2; // Adjust the number of iterations as needed
ArrayList<PVector> vertices;
ArrayList<int[]> faces;

void setup() {
  size(800, 800, P3D);
  PeasyCam cam = new PeasyCam(this, 500);
  initializeMesh();
}

void draw() {
  background(255);
  stroke(0);
  fill(150, 150, 255);
  drawMesh();
}

void initializeMesh() {
  // Initialize your mesh vertices and faces here
  // For example, a cube:
  vertices = new ArrayList<>();
  faces = new ArrayList<>();
  
  //vertices.add(new PVector(0, -50, 0)); // Apex
  //vertices.add(new PVector(-50, 50, -50)); // Base vertex 1
  //vertices.add(new PVector(50, 50, -50));  // Base vertex 2
  //vertices.add(new PVector(50, 50, 50));   // Base vertex 3
  //vertices.add(new PVector(-50, 50, 50));  // Base vertex 4
  
  // Define faces using vertex indices
  //faces.add(new int[]{0, 1, 2}); // Base triangle 1
  //faces.add(new int[]{0, 2, 3}); // Base triangle 2
  //faces.add(new int[]{0, 3, 4}); // Base triangle 3
  //faces.add(new int[]{0, 4, 1}); // Base triangle 4
  //faces.add(new int[]{1, 2, 3, 4}); // Side face
  vertices.add(new PVector(-50, -50, -50));
  vertices.add(new PVector(50, -50, -50));
  vertices.add(new PVector(50, 50, -50));
  vertices.add(new PVector(-50, 50, -50));
  vertices.add(new PVector(-50, -50, 50));
  vertices.add(new PVector(50, -50, 50));
  vertices.add(new PVector(50, 50, 50));
  vertices.add(new PVector(-50, 50, 50));
  
  //// Define faces using vertex indices
  faces.add(new int[]{0, 1, 2, 3});
  faces.add(new int[]{4, 5, 6, 7});
  faces.add(new int[]{0, 3, 7, 4});
  faces.add(new int[]{1, 2, 6, 5});
  faces.add(new int[]{3, 2, 6, 7});
  faces.add(new int[]{0, 1, 5, 4});
}

void drawMesh() {
  beginShape(QUADS);
  for (int[] face : faces) {
    for (int vertexIndex : face) {
      PVector v = vertices.get(vertexIndex);
      vertex(v.x, v.y, v.z);
    }
  }
  endShape();
}

void keyPressed() {
  // Vérifie si la touche "c" est enfoncée
  if (key == 'c' || key == 'C') {
    // Appelle la fonction CatmullClark pour effectuer la subdivision
    catmullClarkSubdivision();
  }
}

void catmullClarkSubdivision() {
  ArrayList<PVector> newVertices = new ArrayList<>();
  ArrayList<int[]> newFaces = new ArrayList<>();

  // Step 1: Calculate face points
  ArrayList<PVector> facePoints = calculate_face_points(vertices, faces);

  // Step 2: Calculate edge points
  HashMap<String, PVector> edgePoints = calculate_edge_points(vertices, faces, facePoints);

  // Step 3: Update vertex coordinates
  for (int i = 0; i < vertices.size(); i++) {
    PVector oldCoords = vertices.get(i);
    ArrayList<PVector> avgFacePoints = getAdjacentFacePoints(i, faces, facePoints);
    ArrayList<PVector> avgMidEdges = getAdjacentEdgeMidpoints(i, vertices, faces, edgePoints);

    int n = avgFacePoints.size();
    float m1 = (n - 3.0) / n;
    float m2 = 1.0 / n;
    float m3 = 2.0 / n;

    PVector newCoords = new PVector();
    newCoords.add(PVector.mult(oldCoords, m1));
    PVector avg = new PVector();
    for (PVector avgFacePoint : avgFacePoints) {
      avg.add(avgFacePoint);
    }
    avg.div(n);
    newCoords.add(PVector.mult(avg, m2));
    avg = new PVector();
    for (PVector avgMidEdge : avgMidEdges) {
      avg.add(avgMidEdge);
    }
    avg.div(n);
    newCoords.add(PVector.mult(avg, m3));
    newVertices.add(newCoords);
  }
  
  // Add face points and edge points to the new vertices
  newVertices.addAll(facePoints);
  for (PVector edgePoint : edgePoints.values()) {
    newVertices.add(edgePoint);
  }
  
  // Step 4: Generate new faces
  for (int[] face : faces) {
    int a = face[0];
    int b = face[1];
    int c = face[2];
    int d = face.length == 4 ? face[3] : -1; // Check if it's a quad or triangle
    int edgePointCD =0;
    int edgePointDA =0;
    int edgePointAB = getEdgePointIndex(a, b, edgePoints,newVertices);
    int edgePointCA = getEdgePointIndex(c, a, edgePoints,newVertices);
    int edgePointBC = getEdgePointIndex(b, c, edgePoints,newVertices);
    if (d != -1) {
      edgePointDA = getEdgePointIndex(d, a, edgePoints,newVertices);
      edgePointCD = getEdgePointIndex(c, d, edgePoints,newVertices);
    }
    

    int facePointIndex = getFacePointIndex(facePoints, face,newVertices);

    int[] newFace1 = new int[4];
    int[] newFace2 = new int[4];;
    int[] newFace3 = new int[4];;
    int[] newFace4 = new int[4];;
    // Create new faces
    if (d == -1) {
      int[] save = {a, edgePointAB, facePointIndex, edgePointCA};
      System.arraycopy(save,0,newFace1,0,save.length);
      int[] save1 = {b, edgePointBC, facePointIndex, edgePointAB};
      System.arraycopy(save1,0,newFace2,0,save1.length);
      int[] save2 = {c, edgePointCA, facePointIndex, edgePointBC};
      System.arraycopy(save2,0,newFace3,0,save2.length);
      //newFace1 = {a, edgePointAB, facePointIndex, edgePointCA};
      //newFace2 = {b, edgePointBC, facePointIndex, edgePointAB};
      //newFace3 = {c, edgePointCA, facePointIndex, edgePointBC};
    }else{
      int[] save3 = {a, edgePointAB, facePointIndex, edgePointDA};
      System.arraycopy(save3,0,newFace1,0,save3.length);
      int[] save4 = {b, edgePointBC, facePointIndex, edgePointAB};
      System.arraycopy(save4,0,newFace2,0,save4.length);
      int[] save5 = {c, edgePointCD, facePointIndex, edgePointBC};
      System.arraycopy(save5,0,newFace3,0,save5.length);
      int[] save6 = {d, edgePointDA, facePointIndex, edgePointCD};
      System.arraycopy(save6,0,newFace4,0,save6.length);
      //newFace1 = {a, edgePointAB, facePointIndex, edgePointDA};
      //newFace2 = {b, edgePointBC, facePointIndex, edgePointAB};
      //newFace3 = {c, edgePointCD, facePointIndex, edgePointBC};
      //newFace4 = {d, edgePointDA, facePointIndex, edgePointCD};
    }

    newFaces.add(newFace1);
    newFaces.add(newFace2);
    newFaces.add(newFace3);

    if (d != -1) {
      newFaces.add(newFace4);
    }
  }

  // Update mesh with new vertices and faces
  vertices = newVertices;
  faces = newFaces;
}

ArrayList<PVector> getAdjacentFacePoints(int vertexIndex, ArrayList<int[]> faces, ArrayList<PVector> facePoints) {
  ArrayList<PVector> adjacentFacePoints = new ArrayList<>();
  for (int[] face : faces) {
    for (int faceVertex : face) {
      if (faceVertex == vertexIndex) {
        adjacentFacePoints.add(facePoints.get(faces.indexOf(face)));
        break;
      }
    }
  }
  return adjacentFacePoints;
}

ArrayList<PVector> getAdjacentEdgeMidpoints(int vertexIndex, ArrayList<PVector> vertices, ArrayList<int[]> faces, HashMap<String, PVector> edgePoints) {
  ArrayList<PVector> adjacentEdgeMidpoints = new ArrayList<>();
  PVector vertex = vertices.get(vertexIndex);

  for (int[] face : faces) {
    for (int i = 0; i < face.length; i++) {
      int nextIndex = (i + 1) % face.length;
      int edgeStart = face[i];
      int edgeEnd = face[nextIndex];

      if ((edgeStart == vertexIndex || edgeEnd == vertexIndex) && edgeStart != edgeEnd) {
        String edgeKey = getEdgeKey(edgeStart, edgeEnd);
        PVector edgeMidpoint = edgePoints.get(edgeKey);
        adjacentEdgeMidpoints.add(edgeMidpoint);
        break;
      }
    }
  }

  return adjacentEdgeMidpoints;
}

int getFacePointIndex(ArrayList<PVector> facePoints, int[] face, ArrayList<PVector> newVertices) {
  PVector avgFacePoint = new PVector();
  for (int vertexIndex : face) {
    avgFacePoint.add(vertices.get(vertexIndex));
  }
  avgFacePoint.div(face.length);
  
  int existingIndex = newVertices.indexOf(avgFacePoint);
  if (existingIndex == -1) {
    newVertices.add(avgFacePoint);
    return newVertices.size() - 1;
  } else {
    return existingIndex;
  }
}


int getEdgePointIndex(int start, int end, HashMap<String, PVector> edgePoints, ArrayList<PVector> newVertices) {
  String edgeKey = getEdgeKey(start, end);
  PVector edgeMidpoint = edgePoints.get(edgeKey);
  
  if (edgeMidpoint == null) {
    edgeMidpoint = new PVector((vertices.get(start).x + vertices.get(end).x) / 2,
                               (vertices.get(start).y + vertices.get(end).y) / 2,
                               (vertices.get(start).z + vertices.get(end).z) / 2);
    edgePoints.put(edgeKey, edgeMidpoint);
  }

  int edgePointIndex = newVertices.indexOf(edgeMidpoint);
  if (edgePointIndex == -1) {
    newVertices.add(edgeMidpoint);
    return newVertices.size() - 1;
  } else {
    return edgePointIndex;
  }
}


String getEdgeKey(int start, int end) {
  return start < end ? start + "-" + end : end + "-" + start;
}

ArrayList<PVector> calculate_face_points(ArrayList<PVector> vertices, ArrayList<int[]> faces) {
  ArrayList<PVector> facePoints = new ArrayList<>();
  for (int[] face : faces) {
    PVector avgPoint = new PVector();
    for (int vertexIndex : face) {
      avgPoint.add(vertices.get(vertexIndex));
    }
    avgPoint.div(face.length);
    facePoints.add(avgPoint);
  }
  return facePoints;
}

HashMap<String, PVector> calculate_edge_points(ArrayList<PVector> vertices, ArrayList<int[]> faces, ArrayList<PVector> facePoints) {
  HashMap<String, PVector> edgePoints = new HashMap<>();
  for (int[] face : faces) {
    for (int i = 0; i < face.length; i++) {
      int nextIndex = (i + 1) % face.length;
      int edgeStart = face[i];
      int edgeEnd = face[nextIndex];

      PVector facePointStart = facePoints.get(faces.indexOf(face));
      PVector facePointEnd = facePoints.get(faces.indexOf(getAdjacentFace(face, vertices, faces, edgeStart, edgeEnd)));

      PVector edgePoint = new PVector();
      edgePoint.add(vertices.get(edgeStart));
      edgePoint.add(vertices.get(edgeEnd));
      edgePoint.add(facePointStart);
      edgePoint.add(facePointEnd);
      edgePoint.div(4.0);

      String edgeKey = getEdgeKey(edgeStart, edgeEnd);
      edgePoints.put(edgeKey, edgePoint);
    }
  }
  return edgePoints;
}

int[] getAdjacentFace(int[] face, ArrayList<PVector> vertices, ArrayList<int[]> faces, int edgeStart, int edgeEnd) {
  for (int[] adjacentFace : faces) {
    if (adjacentFace != face && hasEdge(adjacentFace, edgeStart, edgeEnd)) {
      return adjacentFace;
    }
  }
  return null;
}

boolean hasEdge(int[] face, int edgeStart, int edgeEnd) {
  for (int i = 0; i < face.length; i++) {
    int nextIndex = (i + 1) % face.length;
    if ((face[i] == edgeStart && face[nextIndex] == edgeEnd) ||
        (face[i] == edgeEnd && face[nextIndex] == edgeStart)) {
      return true;
    }
  }
  return false;
}
