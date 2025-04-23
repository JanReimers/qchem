#include <gtkmm.h>
#include "ControllerWindow.H"
#include <iostream>

Glib::RefPtr<Gtk::Application> app;
ControllerWindow* controller;

void on_app_activate()
{
    auto refBuilder = Gtk::Builder::create();
    try
    {
      refBuilder->add_from_file("qchem-gtk.ui");
    }
    catch(const Glib::FileError& ex)
    {
      std::cerr << "FileError: " << ex.what() << std::endl;
      return;
    }
    catch(const Glib::MarkupError& ex)
    {
      std::cerr << "MarkupError: " << ex.what() << std::endl;
      return;
    }
    catch(const Gtk::BuilderError& ex)
    {
      std::cerr << "BuilderError: " << ex.what() << std::endl;
      return;
    }
  
    // Controller c(refBuilder);
    // auto main = refBuilder->get_widget<Gtk::Window>("main");
    controller= Gtk::Builder::get_widget_derived<ControllerWindow>(refBuilder, "main");
    app->add_window(*controller);
    controller->set_visible(true);
}

void on_app_shutdown()
{
  delete controller;
}

int main(int argc, char** argv)
{
    app = Gtk::Application::create("org.gtkmm.example");

    app->signal_activate().connect([] () { on_app_activate(); });
    app->signal_shutdown().connect([] () { on_app_shutdown(); });

  return app->run(argc, argv);
}