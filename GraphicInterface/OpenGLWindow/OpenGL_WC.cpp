#include "OpenGL_WC.h"
#include "WindowConstructor.h"

#include "WeakForm.h"
#include "GmshTools.h"

void Playback::InitializeGL() {
    glClearColor(1.0f, 1.0f, 1.0f, 1); // BackGround Color
    glEnable(GL_DEPTH_TEST); // Enable depth testing. Helps determine which objects are in front of others when they are drawn on the screen.
}

void Playback::draw() {
    static bool firstTime = true;
    if (firstTime) {
        InitializeGL();
        firstTime = false;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clears the color and depth buffer

    // Get current window size
    int w = this->w();
    int h = this->h();

    // Adjust the viewport
    glViewport(0, 0, w, h);

    // view transformations
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = static_cast<float>(w) / static_cast<float>(h);
    gluPerspective(45.0, aspect, 1.0, 1000.0); // Perspective projection with aspect ratio

    //Since gluLookAt defines the center of the camera, lets calculate a average coordenates from geometry:
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, zoom, //Camera Position
              0, 0, 0, //New Camera Center
              0, 1, 0); // Up vector

    //glTranslatef(-average_x, -average_y, 0); //Here , the new camera center.
    
    // Applys rotations
    glRotatef(rotationX, 1, 0, 0);
    glRotatef(rotationY, 0, 1, 0);

    // Draws the geometry
    DrawGeometry(mat);

    //DrawGeometry(mat);
    
    // Draws the axes
    DrawAxes();

}

void Playback::setInt(int newValue) {
    mat = newValue;
    redraw();  // Request a redraw to reflect the new value
}


void Playback::PrintElement(Element *el, CompMesh *cmesh){

    auto connect = el->getConnectivity();
    switch (el->PrintType())
    {
    case 3://Linear line
        {
            auto core_A = cmesh->NodeVec()[connect[0]]->getCoordinates();
            auto core_B = cmesh->NodeVec()[connect[1]]->getCoordinates();

            //Since glVertex3fv requires Glfloats as entry, let's convert to it:
            GLfloat pa []= {static_cast<float>(core_A[0]), static_cast<float>(core_A[1]),static_cast<float> (core_A[2])};
            GLfloat pb []= {static_cast<float>(core_B[0]), static_cast<float>(core_B[1]),static_cast<float> (core_B[2])};
            

            sum_x += core_A[0] + core_B[0];
            sum_y += core_A[1] + core_B[1];
            sum_z += core_A[2] + core_B[2];

            //Creates OpenGL points
            glVertex3fv(pa);
            glVertex3fv(pb);
        }
        break;
        case 5: //Linear triangle
        {
            auto core_A = cmesh->NodeVec()[connect[0]]->getCoordinates();
            auto core_B = cmesh->NodeVec()[connect[1]]->getCoordinates();
            auto core_C = cmesh->NodeVec()[connect[2]]->getCoordinates();

            //Since glVertex3fv requires Glfloats as entry, let's convert to it:
            GLfloat pa []= {static_cast<float>(core_A[0]), static_cast<float>(core_A[1]),static_cast<float> (core_A[2])};
            GLfloat pb []= {static_cast<float>(core_B[0]), static_cast<float>(core_B[1]),static_cast<float> (core_B[2])};
            GLfloat pc []= {static_cast<float>(core_C[0]), static_cast<float>(core_C[1]),static_cast<float> (core_C[2])};

            //Creates OpenGL points
            glVertex3fv(pa);
            glVertex3fv(pb);

            glVertex3fv(pb);
            glVertex3fv(pc);

            glVertex3fv(pc);
            glVertex3fv(pa);
        }
        break;
        
        case 9: //linear quadrilateral
        {
            auto core_A = cmesh->NodeVec()[connect[0]]->getCoordinates();
            auto core_B = cmesh->NodeVec()[connect[1]]->getCoordinates();
            auto core_C = cmesh->NodeVec()[connect[2]]->getCoordinates();
            auto core_D = cmesh->NodeVec()[connect[3]]->getCoordinates();

            //Since glVertex3fv requires Glfloats as entry, let's convert to it:
            GLfloat pa []= {static_cast<float>(core_A[0]), static_cast<float>(core_A[1]),static_cast<float> (core_A[2])};
            GLfloat pb []= {static_cast<float>(core_B[0]), static_cast<float>(core_B[1]),static_cast<float> (core_B[2])};
            GLfloat pc []= {static_cast<float>(core_C[0]), static_cast<float>(core_C[1]),static_cast<float> (core_C[2])};
            GLfloat pd []= {static_cast<float>(core_D[0]), static_cast<float>(core_D[1]),static_cast<float> (core_D[2])};

            //Creates OpenGL points
            glVertex3fv(pa);
            glVertex3fv(pb);

            glVertex3fv(pb);
            glVertex3fv(pc);

            glVertex3fv(pc);
            glVertex3fv(pd);

            glVertex3fv(pd);
            glVertex3fv(pa);
        }
        break;

        case 10: //linear tetrahedron
        {
            auto core_A = cmesh->NodeVec()[connect[0]]->getCoordinates();
            auto core_B = cmesh->NodeVec()[connect[1]]->getCoordinates();
            auto core_C = cmesh->NodeVec()[connect[2]]->getCoordinates();
            auto core_D = cmesh->NodeVec()[connect[3]]->getCoordinates();

            //Since glVertex3fv requires Glfloats as entry, let's convert to it:
            GLfloat pa []= {static_cast<float>(core_A[0]), static_cast<float>(core_A[1]),static_cast<float> (core_A[2])};
            GLfloat pb []= {static_cast<float>(core_B[0]), static_cast<float>(core_B[1]),static_cast<float> (core_B[2])};
            GLfloat pc []= {static_cast<float>(core_C[0]), static_cast<float>(core_C[1]),static_cast<float> (core_C[2])};
            GLfloat pd []= {static_cast<float>(core_D[0]), static_cast<float>(core_D[1]),static_cast<float> (core_D[2])};

            //Creates OpenGL points
            glVertex3fv(pa);
            glVertex3fv(pb);

            glVertex3fv(pb);
            glVertex3fv(pc);

            glVertex3fv(pc);
            glVertex3fv(pd);

            glVertex3fv(pd);
            glVertex3fv(pa);

            
        }
        break;


    default:
        PanicButton();
        break;
    }
    
}



void Playback::DrawGeometry(int matidcolor) {
    if(matids.size() == 0){
        return;
    }

    CompMesh* Outmesh = new CompMesh();
    for(int i = 0; i < matids.size(); i++){
        WeakForm wf;
        wf.Id() = matids[i];
        Outmesh->InsertMaterial(&wf);
    }

    GmshTools::Read(*Outmesh, mshpath);
    glBegin(GL_LINES);
    //Sets line thickness
    glLineWidth(2.0f);

    //Renders the Geometry (lines)
    sum_x = 0.;
    sum_y = 0.;
    sum_z = 0.;

    for (int i = 0; i < Outmesh->NElements(); i++)
    {
        auto geo = Outmesh->ElementVec()[i];
        
        if (geo->GetWeakForm()->Id() == matidcolor){
            
            glColor3f(1,0,0);
            PrintElement(geo,Outmesh);

        } else {

            glColor3f(0,0,0);
            PrintElement(geo,Outmesh);
            
        }
    }
    average_x = sum_x/Outmesh->NNodes();
    average_y = sum_y/Outmesh->NNodes();
    average_z = sum_z/Outmesh->NNodes();
    delete Outmesh;

    glEnd();
}


void Playback::DrawAxes() {
     // Set viewport for the axes in the corner
    glViewport(10, 10, 100, 100);

    // Set up orthogonal projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);

    // Set modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Apply the same rotations as the main view
    glRotatef(rotationX, 1, 0, 0);
    glRotatef(rotationY, 0, 1, 0);

    // Set line width for the axes within the supported range
    glLineWidth(2.0f); // Adjust this value within the supported range

    // Draw the axes
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);

    // X axis
    glVertex3f(0, 0, 0);
    glVertex3f(0.5, 0, 0);

    // Y axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0.5, 0);

    // Z axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 0.5);
    glEnd();

    // Draw arrows at the end of each axis
    glPushMatrix();
    glTranslatef(0.5, 0, 0);
    glutSolidCone(0.05, 0.1, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, 0.5, 0);
    glutSolidCone(0.05, 0.1, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, 0, 0.5);
    glutSolidCone(0.05, 0.1, 10, 10);
    glPopMatrix();

    // Draw labels for each axis
    glRasterPos3f(0.55, 0, 0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'X');

    glRasterPos3f(0, 0.55, 0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Y');

    glRasterPos3f(0, 0, 0.55);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'Z');

    // Draw a 3D sphere at the origin
    // Set viewport for the sphere in the center (0,0,0)
    glViewport(0, 0, 714, 703);
    glPushMatrix();
    //glTranslatef(0.0f, 0.0f, 0.0f); // Explicitly set the position to the origin
    glColor3f(0, 0, 1); // Set the sphere color (RGB)
    glutSolidSphere(0.01, 50, 50); // Draw a sphere with radius 0.1
    glPopMatrix();
}

int Playback::handle(int event) {
    switch (event) {
    case FL_PUSH:
        last_x = Fl::event_x();
        last_y = Fl::event_y();
        rotating = true;
        return 1;

    case FL_DRAG:
        if (rotating) {
            int dx = Fl::event_x() - last_x;
            int dy = Fl::event_y() - last_y;

            rotationX += dy; // Update rotation angle for X axis
            rotationY += dx; // Update rotation angle for Y axis

            last_x = Fl::event_x();
            last_y = Fl::event_y();

            redraw();
            return 1;
        }
        break;

    case FL_RELEASE:
        rotating = false;
        return 1;

    case FL_MOUSEWHEEL:
        zoom += Fl::event_dy() * 5.f; // Adjust zoom based on scroll direction
        redraw();
        return 1;

    default:
        break;
    }

    return Fl_Gl_Window::handle(event);
}

void IdleCallback(void* pData) {
    if (pData != NULL) {
        Playback* pWindow = reinterpret_cast<Playback*>(pData);
        if (pWindow->animating) {
            pWindow->rotation += pWindow->rotationIncrement / 100;
            pWindow->redraw();
        }
    }
}

Playback::Playback(int X, int Y, int W, int H, const char *L) 
    : Fl_Gl_Window(X, Y, W, H, L), last_x(0), last_y(0), rotating(false), 
      rotationX(0), rotationY(0), zoom(8) { // Initialize zoom to -1.0 to start further back
    mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE);
    rotation = 0.f;
    rotationIncrement = 10.f;
    animating = false;
    Fl::add_idle(IdleCallback, this);
    end();
}
