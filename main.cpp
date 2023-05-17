// 2017-09-27 testing for potential for another bmk

#include <vector>
#include <unordered_set>
#include <map>
#include <glm/glm.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "bvh.hpp"
#include "utils.hpp"

#ifdef VISUALIZE
	#include <stdlib.h>
	#include <stdio.h>
	#include <GL/glew.h>
	#ifdef __APPLE__
		#include <GLUT/glut.h>
		#include <OpenGL/glu.h>
	#else
		#include <GL/freeglut.h>
		#include <GL/glu.h>
	#endif
	#include "openglstuff.hpp"
#endif

// CUDA
BVH BVHGen_GPU(const std::vector<Vec3>&, const std::vector<Vec3>&);
bool g_use_cuda = false;
int  g_bvh_method = 1;

#ifdef USE_CUDA
extern int g_num_blk, g_num_thd;
#endif

int main(int argc, char** argv) {
	{
		char* x;
		x = getenv("TESTFLAG");
		if (x) g_bvh_method = 2;
		x = getenv("USE_CUDA");
		if (x) g_use_cuda = (bool)atoi(x);
	}
	char filename[100];
	sprintf(filename, "dragon.obj");
	for (int i=1; i<argc; i++) {
		#ifdef USE_CUDA
		if (2 == sscanf(argv[i], "dim=%d,%d", &g_num_blk, &g_num_thd))
			printf("Dim changed to <<<%d, %d>>>\n", g_num_blk, g_num_thd);
		#endif
		if (2 == sscanf(argv[i], "model=%s", filename))
			printf("Model changed to %s\n", filename);
	}

	std::vector<Vec3> verts, centroids_transformed;
	std::vector<std::vector<int> > idxes;
	FILE* f = fopen(filename, "r");
	if (!f) { printf("ERROR: File %s does not exist\n", filename); assert(0); }
	char* line = (char*)malloc(233); size_t nbytes;
	int line_idx = 0;
	while (int read = getline(&line, &nbytes, f) > 0) {
		if (line[0] == 'v') {
			Vec3 v;
			if (not (3 == sscanf(line, "v %f %f %f", &v.x, &v.y, &v.z))) {
				printf("Line %d\n", line_idx + 1);
				assert (0);
			}
			verts.push_back(v);
		} else if (line[0] == 'f') {
			int i0, i1, i2;
			if (not (3 == sscanf(line, "f %d %d %d", &i0, &i1, &i2))) {
				printf("Line %d\n", line_idx + 1);
				assert (0);
			}
			idxes.push_back({i0, i1, i2});
		} else { }
		line_idx ++;
	};
	printf("%lu vertices, %lu faces\n", verts.size(), idxes.size());

	Vec3 c_min(1e20, 1e20, 1e20), c_max(-1e20, -1e20, -1e20);
	for (int i=0; i<int(idxes.size()); i++) {
		Vec3 p0 = verts[idxes[i][0]-1], p1 = verts[idxes[i][1]-1],
		          p2 = verts[idxes[i][2]-1],
				  centroid = (p0 + p1 + p2) * (1.0f / 3.0f);
		c_min = Vec3Min(centroid, c_min);
		c_max = Vec3Max(centroid, c_max);
		centroids_transformed.push_back(centroid);
	}

	std::vector<Vec3> centroids = centroids_transformed;

	for (int i=0; i<int(centroids_transformed.size()); i++) {
		Vec3* p = &(centroids_transformed[i]);
		p->x = (p->x - c_min.x) / (c_max.x - c_min.x);
		p->y = (p->y - c_min.y) / (c_max.y - c_min.y);
		p->z = (p->z - c_min.z) / (c_max.z - c_min.z);
	}

	BVH bvh(centroids);
	bvh.points_transformed = centroids_transformed;
	if (g_bvh_method == 1) bvh.Build();
	else bvh.Build_Karras();

	#ifdef USE_CUDA
	if (g_use_cuda) { 
		bvh = BVHGen_GPU(centroids, centroids_transformed);
	}
	#endif

	#ifdef VISUALIZE
	const int W = 800, H = 480;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH);
	#ifndef __APPLE__
	glutInitContextVersion(1, 3);
	#endif
	glutInitWindowSize(W, H);
	glutInit(&argc, argv);
	glutCreateWindow("BVH Construction");
	glutDisplayFunc(render);
	glutIdleFunc(update);
	glutKeyboardFunc(keydown);
	glutKeyboardUpFunc(keyup);
	glutSpecialFunc(keydown2);
	if(glewInit() == GLEW_OK) {

		// The model itself
		std::vector<std::vector<int> > mapped_idxes;
		for (int i=0; i<bvh.mc.size(); i++) {
			const unsigned code = bvh.mc.at(i);
			for (const int ii : bvh.mc2point.at(code)) {
				mapped_idxes.push_back(idxes.at(ii));
			}
		}
		g_mesh = new Mesh(verts, mapped_idxes);
		g_bvh_vis = new BVH_Vis(bvh);

		// The BVH

		// Initialize viewport
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(60, (float)W/(float)H, 0.1, 500);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glutMainLoop();
	}
	#endif

	return 0;
}
