package ballista;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;
import java.util.Set;

/**
 *
 * @author tsunemat
 */

public class ReadInit {

    private final Properties prop;

    public ReadInit(Properties prop) {
        this.prop = prop;
    }

    public synchronized int size() {
        return prop.size();
    }

    public synchronized boolean containsKey(String key) {
        return prop.containsKey(key);
    }

    public Set<String> keySet() {
        return prop.stringPropertyNames();
    }

    public String getString(String key) {
        return prop.getProperty(key);
    }

    public int getInt(String key) {
        return Integer.parseInt(prop.getProperty(key));
    }

    public double getDouble(String key) {
        return Double.parseDouble(prop.getProperty(key));
    }

    public boolean getBoolean(String key) {
        return Boolean.parseBoolean(prop.getProperty(key));
    }


}

