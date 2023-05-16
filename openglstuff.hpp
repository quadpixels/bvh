void MyCheckGLError(const char* tag) {
	GLenum err = glGetError();
	if (err != GL_NO_ERROR) {
		printf("[%s] OpenGL Error: %d\n", tag, err);
		assert(0);
	}
}

class Mesh {
public:
	unsigned vertex_vbo, index_vbo, index_size;
	Mesh(const std::vector<Vec3>& verts,
	     const std::vector<std::vector<int> >& faces);
	void Render();
	void RenderRange(int idx_lb, int idx_ub);
};

Mesh* g_mesh = NULL;

Mesh::Mesh(const std::vector<Vec3>& verts,
	     const std::vector<std::vector<int> >& faces) {
	printf("Mesh::Mesh called |verts|=%lu |faces|=%lu\n",
		verts.size(), faces.size());
	unsigned buffers[2];
	glGenBuffers(2, buffers);
	MyCheckGLError("Gen Buffer");

	int idx = 0;

	vertex_vbo = buffers[0];
	index_vbo  = buffers[1];

	float* tmpf = new float[verts.size() * 3];
	for (int i=0; i<verts.size(); i++) {
		for (int j=0; j<3; j++) {
			tmpf[idx] = verts[i][j];
			idx ++;
	} }
	glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo);
	glBufferData(GL_ARRAY_BUFFER, idx*sizeof(float), tmpf, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	MyCheckGLError("Vertex VBO");
	
	idx = 0;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
	int* tmp = new int[faces.size() * 3];
	for (int i=0; i<faces.size(); i++) {
		for (int j=0; j<3; j++) {
			tmp[idx] = faces[i][j] - 1; // MUST MINUS ONE
			idx ++;
		}
	}
	index_size = idx;
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx*sizeof(int), tmp, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MyCheckGLError("Index VBO");
	delete[] tmp;
	delete[] tmpf;
}

void Mesh::Render() { RenderRange(0, index_size-1); }

void Mesh::RenderRange(int idx_lb, int idx_ub) {
	glEnableClientState(GL_VERTEX_ARRAY);
	glColor3f(0.7f, 0.7f, 0.7f);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
	const int s = idx_ub - idx_lb + 1;
	int arg0, arg1;
	GLvoid* arg2;
	arg0 = idx_lb;
	arg1 = idx_ub-1;
	arg2 = (GLvoid*)(sizeof(int)*idx_lb);

	// Range 指的是值的 range（意同range query）,所以不能用

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, s, GL_UNSIGNED_INT, arg2);
	MyCheckGLError("draw elements");

	glPolygonMode(GL_FRONT, GL_LINE);
	glColor3f(0.3f, 0.3f, 0.3f);
	glDrawElements(GL_TRIANGLES, s, GL_UNSIGNED_INT, arg2);

	glPolygonMode(GL_BACK, GL_LINE);
	glColor3f(0.3f, 0.0f, 0.0f);
	glDrawElements(GL_TRIANGLES, s, GL_UNSIGNED_INT, arg2);

	glBindBuffer(GL_ARRAY_BUFFER,         0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

class BVH_Vis {
public:
	BVH* bvh;
	int vertex_vbo, index_vbo, curr_node_idx;
	BVH_Vis(BVH& bvh);
	void Render();
	int GetNeighbor(char dir);
	void MoveToLeft();
	void MoveToRight();
	void MoveToParent();
	void PrintSummary();
};

BVH_Vis::BVH_Vis(BVH& bvh) {
	this->bvh = &bvh;
	curr_node_idx = 0;
	int nc = int(bvh.nodes.size());
	unsigned buffers[2];
	glGenBuffers(2, buffers);
	
	vertex_vbo = buffers[0];
	index_vbo  = buffers[1];

	float* v   = new float[nc*3 * 8]; // 8 Vertices per BVH Node
	int* idxes = new int[nc*24];     // 12 indices per BVH Node
	for (int i=0; i<nc; i++) {
		Vec3 p0 = bvh.nodes[i].min_ext,
		          p1 = bvh.nodes[i].max_ext;
		float x0 = p0.x, y0 = p0.y, z0 = p0.z,
		      x1 = p1.x, y1 = p1.y, z1 = p1.z;
		float aabb[] = {
			x0, y0, z0,   x0, y0, z1,   x1, y0, z1,   x1, y0, z0,
			x0, y1, z0,   x0, y1, z1,   x1, y1, z1,   x1, y1, z0,
		};
		int aabbidx[] = { 
			0,1,1,2,2,3,3,0,
			0,4,1,5,2,6,3,7,
			4,5,5,6,6,7,7,4,
		};
		for (int j=0; j<24; j++) v[i*24+j]     = aabb[j];
		for (int j=0; j<24; j++) idxes[i*24+j] = aabbidx[j] + 8*i;
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
	MyCheckGLError("111");
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, nc*24*sizeof(int), idxes, GL_STATIC_DRAW);
	MyCheckGLError("112");
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MyCheckGLError("113");

	glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo);
	glBufferData(GL_ARRAY_BUFFER, nc*24*sizeof(float), v, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	MyCheckGLError("Gen VBO for BVH Visualization");
	delete[] v;
	delete[] idxes;
}

void BVH_Vis::Render() {
	glPolygonMode(GL_FRONT, GL_LINE);
	glDisable(GL_LIGHTING);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glColor3f(0.0f, 0.5f, 0.0f);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
	int lb = 24*curr_node_idx, ub = lb + 23;
	glDrawRangeElements(GL_LINES, 0, 23, 24, GL_UNSIGNED_INT,
		(GLvoid*)(sizeof(int)*24*curr_node_idx));

	glBindBuffer(GL_ARRAY_BUFFER,         0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glPolygonMode(GL_FRONT, GL_FILL);

	MyCheckGLError("draw bounding volume");
}

int BVH_Vis::GetNeighbor(char dir) {
	BVH::BVHNode* curr_node = &(bvh->nodes[curr_node_idx]);
	switch (dir) {
		case 'L': return curr_node->left;
		case 'R': return curr_node->right;
		case 'P': return curr_node->parent;
	}
	return -999;
}

void BVH_Vis::MoveToLeft() {
	int x = GetNeighbor('L');
	if (x != -999) curr_node_idx = x;
}

void BVH_Vis::MoveToRight() {
	int x = GetNeighbor('R');
	if (x != -999) curr_node_idx = x;
}

void BVH_Vis::MoveToParent() {
	int x = GetNeighbor('P');
	if (x != -999) curr_node_idx = x;
}

void BVH_Vis::PrintSummary() {
	BVH::BVHNode* curr_node = &(bvh->nodes[curr_node_idx]);
	printf("Now showing BVHNode[%d]=%d faces [%d,%d]=(%g,%g,%g)-->(%g,%g,%g)\n",
		curr_node_idx,
		curr_node->node_ub - curr_node->node_lb + 1,
		curr_node->node_lb,
		curr_node->node_ub,
		curr_node->min_ext.x, curr_node->min_ext.y, curr_node->min_ext.z,
		curr_node->max_ext.x, curr_node->max_ext.y, curr_node->max_ext.z
	);
}

BVH_Vis* g_bvh_vis = NULL;

Vec3 g_crysball(0, -5, -15), g_crysball_v(0, 0, 0);

void keydown(unsigned char key, int x, int y) {
	switch (key) {
		case 'w': g_crysball_v.z =  0.1f; break;
		case 's': g_crysball_v.z = -0.1f; break;
		case 'a': g_crysball_v.x =  0.1f; break;
		case 'd': g_crysball_v.x = -0.1f; break;
	}
}

void keydown2(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_UP:    
			g_bvh_vis->MoveToParent(); 
			g_bvh_vis->PrintSummary();
		break;
		case GLUT_KEY_LEFT:  
			g_bvh_vis->MoveToLeft(); 
			g_bvh_vis->PrintSummary();
		break;
		case GLUT_KEY_RIGHT: 
			g_bvh_vis->MoveToRight(); 
			g_bvh_vis->PrintSummary();
		break;
	}
}

void keyup(unsigned char key, int x, int y) {
	switch (key) {
		case 'w': 
		case 's': g_crysball_v.z =  0.0f; break;
		case 'a': 
		case 'd': g_crysball_v.x =  0.0f; break;
	}
}

void render() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f, 1.0f, 0.8f, 1.0f);

	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);


	glTranslatef(g_crysball.x, g_crysball.y, g_crysball.z);

	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(5, 0, 0);
	glEnd();

	glColor3f(0, 1, 0);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 5, 0);
	glEnd();

	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 5);
	glEnd();

	if (g_bvh_vis) {
		g_bvh_vis->Render();
		int x = g_bvh_vis->curr_node_idx,
		    lb = g_bvh_vis->bvh->nodes[x].node_lb, // Face IDX
			ub = g_bvh_vis->bvh->nodes[x].node_ub; // Face IDX
		g_mesh->RenderRange(lb*3, ub*3); // Vertex IDX
	} else g_mesh->Render();


	MyCheckGLError("Rendering mesh");

	glutSwapBuffers();
}

void update() {
	g_crysball = g_crysball + g_crysball_v;
	glutPostRedisplay();
}

