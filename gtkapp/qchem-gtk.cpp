#include <gtkmm.h>
#include <iostream>

Glib::RefPtr<Gtk::Application> app;

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
  
    auto main = refBuilder->get_widget<Gtk::Window>("main");
    app->add_window(*main);
    main->set_visible(true);
}

int main(int argc, char** argv)
{
    app = Gtk::Application::create("org.gtkmm.example");

    app->signal_activate().connect([] () { on_app_activate(); });

  return app->run(argc, argv);
}